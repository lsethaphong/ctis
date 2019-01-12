!======================================================================
! FILENAME: sps_SHG.f
!
! DESCRIPTION:
!    This subroutine gives the initial g_k forward projection utilizing
!    the generalized shift invariant basis matrix H
!
! HISTORY:
!   February 7, 2005: Original Code
!
!   Author:           Latsavongsakda Sethaphong
!   Institution:      Molecular Physiology & Biophysics
!                     Piston Lab
!                     Vanderbilt University
!                     Nashville, TN
!                     latsavongsakda.sethapong@vanderbilt.edu
!----------------------------------------------------------------------
       subroutine sps_SHG(shgom,hbasis,g_o,g_k) 

       REAL, DIMENSION(:), INTENT(OUT) :: shgom
       REAL, DIMENSION(:,2), INTENT(IN) :: hbasis
       REAL, DIMENSION(:), INTENT(IN) :: g_o
       REAL, DIMENSION(:), INTENT(OUT) :: g_k
         
       integer, parameter :: r8 = SELECTED_REAL_KIND(15,307)
       integer, parameter :: r4 = SELECTED_REAL_KIND(6,37)
       integer, parameter :: i1 = SELECTED_INT_KIND(2)
       integer, parameter :: i2 = SELECTED_INT_KIND(4)
       integer, parameter :: i4 = SELECTED_INT_KIND(8)
 
       integer, parameter :: MASTER = 0
       integer, parameter :: FROM_MASTER = 1
       integer, parameter :: FROM_WORKER = 2

       include 'mpif.h'
       integer status(MPI_STATUS_SIZE), ierr
       integer numtasks,taskid,numworkers,source,dest,nbytes,mtype
       character*25 g_obs_file, hpinv_nm_file
       integer l-n  
       real(r4), dimension(h_m):: g_p
       real(r4), dimension(:,2):: hshift
       integer, dimension(vxl_cnt,2) :: hrow
!======================================================================
! f_k is shared amongst all, g_obs is broadcast to all
! preallocated
!======================================================================

       call MPI_INIT(ierr)
       call MPI_COMM_RANK(MPI_COMM_WORLD, taskid, ierr)
       call MPI_COMM_SIZE(MPI_COMM_WORLD, numtasks, ierr)
     
       numworkers = numtasks-1
      
C *************************** master task *****************************
       if (taskid .eq. MASTER) then

!======================================================================
!     Send work schedules to remote tasks 
!======================================================================
       mtype = FROM_MASTER
       do 50 j = 1,numworkers
       call MPI_SEND(f_k, h_n, MPI_DOUBLE_PRECISION,
     .    dest, mtype, MPI_COMM_WORLD, ierr)
       call MPI_SEND(h_mn, SIZE(h_mn), MPI_DOUBLE_PRECISION,
     .    dest, mtype, MPI_COMM_WORLD, ierr)
! send section information
       call MPI_SEND(srt_pxl, 1, MPI_INTEGER, dest, mtype,
     .    MPI_COMM_WORLD, ierr)
       call MPI_SEND(stp_pxl, 1, MPI_INTEGER, dest, mtype,
     .    MPI_COMM_WORLD, ierr)

50     continue
!======================================================================
!     Receive g_k results from worker tasks
!     g_k sections
!======================================================================
       mtype = FROM_WORKER
       do 60 i=1, numworkers
       source = i
       call MPI_RECV(g_p, h_m, MPI_DOUBLE_PRECISION,
     .    source, mtype, MPI_COMM_WORLD, status, ierr)
!======================================================================
!gather elements of g_k, just super impose the returned g_p   
!======================================================================
         g_k = g_k + g_p
60    continue

       endif
!======================================================================
!
!                            worker task 
!
!======================================================================
       if (taskid > MASTER) then
C     Receive matrix data from master task
       mtype = FROM_MASTER
       call MPI_RECV(f_k, h_n, MPI_DOUBLE_PRECISION, 
     .       MASTER, mtype, MPI_COMM_WORLD, status, ierr)
       call MPI_RECV(h_mn, SIZE(h_mn), MPI_DOUBLE_PRECISION,
     .       MASTER, mtype, MPI_COMM_WORLD, status, ierr)
       call MPI_RECV(srt_pxl, 1, MPI_INTEGER, MASTER,
     .     mtype, MPI_COMM_WORLD, status, ierr)
       call MPI_RECV(stp_pxl, 1, MPI_INTEGER, 
     .     MASTER, mtype, MPI_COMM_WORLD, status, ierr)

!======================================================================
!
!     G = H+.F               Do matrix multiply 
!
!======================================================================
! 
! need to get pixel id and corresponding row id of f_p
        n = (srt_pxl-1)*vxl_cnt +1
        do 100 k = srt_pxl, stp_pxl ! which pixels
           x_pos = MOD (k,fpa_x)
           y_pos = (k-x_pos)/fpa_x + 1 ! get y position
           if x_pos.eq.0 then ! final check on x position
           x_pos = fpa_x
           endif
           ! is interpreted as row sparse column successive
           call spsMAT_shift(h_mn,x_pos,y_pos,fpa_x,hshift) 
           do 120 m = 1, vxl ! decompse to voxels
              call strip(hrow,hshift,m) ! gets shifted row/col
              n = n + 1
              do 110 l = hrow(m,1),hrow(m,2) ! indices of hpinv
              g_p(INT(abs(hshift(l,2)))  = g_p(INT(abs(hshift(1,2)))  
     ~        + f_k(INT(abs(hshift(l,2)))*hshift(l,1)
110           continue
120        continue 
100     continue

!======================================================================
!   Send basic construct back to master task
!======================================================================

       mtype = FROM_WORKER
       call MPI_SEND(g_p, h_m, MPI_DOUBLE_PRECISION,
     .     MASTER, mtype, MPI_COMM_WORLD, ierr)
       endif
       call MPI_FINALIZE(ierr)

       end program
