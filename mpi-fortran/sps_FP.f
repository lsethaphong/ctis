!======================================================================
! FILENAME: sps_FP.f
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
       subroutine sps_FP(
     &            g_k, ! returns g_k
     &            f_k, ! given f_k
     &            h_mn, ! basis transfer matrix
     &            h_m,h_n,vxl_cnt,fpa_x)

       use com_methods

       INTEGER, INTENT(IN) :: h_m
       INTEGER, INTENT(IN) :: h_n
       INTEGER, INTENT(IN) :: vxl_cnt
       INTEGER, INTENT(IN) :: fpa_x
       REAL, DIMENSION(:), INTENT(IN) :: f_k
       REAL, DIMENSION(:), INTENT(OUT) :: g_k
       REAL, DIMENSION(:,:), INTENT(IN) :: h_mn
         
       integer, parameter :: MASTER = 0
       integer, parameter :: FROM_MASTER = 1
       integer, parameter :: FROM_WORKER = 2

       include 'mpif.h'
       integer status(MPI_STATUS_SIZE), ierr,i,j,k,l,m,n
       integer numtasks,taskid,numworkers,source,dest,nbytes,mtype
       character*25 g_obs_file, hpinv_nm_file
       integer idx, srt_pxl, stp_pxl, extra,avg_pxl
       integer x_pos, y_pos  
       real(r4), dimension(h_m) :: g_p
       real, dimension(SIZE(h_mn,1),SIZE(h_mn,2)) :: hshift
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
       extra = mod(numworkers,h_n/vxl_cnt)
       avg_pxl = (h_n/vxl_cnt-extra)/numworkers
       g_k = 0.
       g_p = 0.
!======================================================================
!     Send work schedules to remote tasks 
!======================================================================
       mtype = FROM_MASTER
       do 50 j = 1,numworkers
          if (j.eq.1) then
             srt_pxl = 1
             stp_pxl = j*avg_pxl + extra
          else
             srt_pxl = (j-1)*avg_pxl + 1 + extra
             stp_pxl = j*avg_pxl + extra
          endif
       call MPI_SEND(f_k, h_n, MPI_DOUBLE_PRECISION,
     .    dest, mtype, MPI_COMM_WORLD, ierr)
       call MPI_SEND(h_mn, SIZE(h_mn), MPI_DOUBLE_PRECISION,
     .    dest, mtype, MPI_COMM_WORLD, ierr)
! send section information
       call MPI_SEND(srt_pxl, 1, MPI_INTEGER, dest, mtype,
     .    MPI_COMM_WORLD, ierr)
       call MPI_SEND(stp_pxl, 1, MPI_INTEGER, dest, mtype,
     .    MPI_COMM_WORLD, ierr)
       call MPI_SEND(fpa_x, 1, MPI_INTEGER, dest, mtype,
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
       call MPI_RECV(fpa_x, 1, MPI_INTEGER,
     .       MASTER, mtype, MPI_COMM_WORLD, status, ierr)

!======================================================================
!
!     G = H+.F               Do matrix multiply 
!
!======================================================================
!      
! need to get pixel id and corresponding row id of f_p
        g_p = 0.
        n = (srt_pxl-1)*vxl_cnt 
        do 100 k = srt_pxl, stp_pxl ! which pixels
           x_pos = MOD(k,fpa_x)
           y_pos = (k-x_pos)/fpa_x + 1 ! get y position
           if (x_pos.eq.0) then ! final check on x position
           x_pos = fpa_x
           endif
           ! is interpreted as row sparse column successive
           call spsMAT_shift(h_mn,x_pos,y_pos,fpa_x,hshift) 
           do 120 m = 1, vxl_cnt ! decompse to voxels
              call strip(hrow,hshift,m) ! gets shifted row/col
              n = n + 1
              do 110 l = hrow(m,1),hrow(m,2) ! indices of hpinv
              idx = INT(abs(hshift(l,2)))
              g_p(idx) = (g_p(idx)+ f_k(idx)*hshift(l,1))
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

       end subroutine sps_FP
