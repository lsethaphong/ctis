!======================================================================
! FILENAME: sps_MAXF.f
!
! DESCRIPTION:
!    
!    
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
       subroutine sps_MAXF(
     ~            f_k, ! f_k is replaced with f_k + 1
     ~            h_mod, ! basis column modulus, indexed by voxel
     ~            hbasis, ! basis forward matrix 
     ~            g_o,  ! g observed
     ~            g_k,  ! g_k back projected
     ~            fpa_x,vxl_cnt)

       use com_methods
       REAL, DIMENSION(:), INTENT(INOUT) :: f_k
       REAL, DIMENSION(:), INTENT(IN) :: g_o
       REAL, DIMENSION(:), INTENT(IN) :: g_k
       REAL, DIMENSION(:), INTENT(IN) :: h_mod
       REAL, DIMENSION(:,:), INTENT(IN) :: hbasis
       INTEGER, INTENT(IN) :: fpa_x
       INTEGER, INTENT(IN) :: vxl_cnt
 
       integer, parameter :: MASTER = 0
       integer, parameter :: FROM_MASTER = 1
       integer, parameter :: FROM_WORKER = 2

       include 'mpif.h'
       integer status(MPI_STATUS_SIZE), ierr
       integer numtasks,taskid,numworkers,source,dest,nbytes,mtype
       integer srt_pxl,stp_pxl,extra,avg_pxl,yi,ye,nf 
       integer i,j,k,l,m,n, x_pos, y_pos
       real tmpsum
       real, dimension(SIZE(g_o,1)):: g_p
       real, dimension(SIZE(hbasis,1),SIZE(hbasis,2)):: hshift
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
       extra = mod(numworkers,INT(h_n/vxl_cnt))
       avg_pxl = (h_n/vxl_cnt-extra)/numworkers
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
        yi =(srt_pxl-1)*vxl_cnt +1
        ye =stp_pxl*vxl_cnt 
       call MPI_SEND(srt_pxl, 1, MPI_INTEGER, dest, mtype,
     .    MPI_COMM_WORLD, ierr)
       call MPI_SEND(stp_pxl, 1, MPI_INTEGER, dest, mtype,
     .    MPI_COMM_WORLD, ierr)
       call MPI_SEND(f_k(yi:ye), (ye-yi+1), MPI_DOUBLE_PRECISION,
     .    dest, mtype, MPI_COMM_WORLD, ierr)
       call MPI_SEND(hbasis, SIZE(hbasis), MPI_DOUBLE_PRECISION,
     .    dest, mtype, MPI_COMM_WORLD, ierr)
       call MPI_SEND(h_mod, vxl_cnt, MPI_DOUBLE_PRECISION,
     .    dest, mtype, MPI_COMM_WORLD, ierr)

50     continue
!======================================================================
!     Receive g_k results from worker tasks
!     g_k sections
!======================================================================
       mtype = FROM_WORKER
       do 60 i=1, numworkers
       source = i
          if (i.eq.1) then
             srt_pxl = 1
             stp_pxl = i*avg_pxl + extra
          else
             srt_pxl = (i-1)*avg_pxl + 1 + extra
             stp_pxl = i*avg_pxl + extra
          endif
        yi =(srt_pxl-1)*vxl_cnt +1
        ye =stp_pxl*vxl_cnt
       call MPI_RECV(f_k(yi:ye), (ye-yi+1), MPI_DOUBLE_PRECISION,
     .    source, mtype, MPI_COMM_WORLD, status, ierr)
!======================================================================
!gather elements of g_k, just super impose the returned g_p   
!======================================================================
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
       call MPI_RECV(srt_pxl, 1, MPI_INTEGER, MASTER,
     .     mtype, MPI_COMM_WORLD, status, ierr)
       call MPI_RECV(stp_pxl, 1, MPI_INTEGER,
     .     MASTER, mtype, MPI_COMM_WORLD, status, ierr)
       call MPI_RECV(f_k, h_n, MPI_DOUBLE_PRECISION, 
     .       MASTER, mtype, MPI_COMM_WORLD, status, ierr)
       call MPI_RECV(hbasis, SIZE(hbasis), MPI_DOUBLE_PRECISION,
     .       MASTER, mtype, MPI_COMM_WORLD, status, ierr)
       call MPI_RECV(h_mod, vxl_cnt, MPI_DOUBLE_PRECISION,
     .       MASTER, mtype, MPI_COMM_WORLD, status, ierr)
!======================================================================
!
!     G = H+.F               Do matrix multiply 
!
!======================================================================
! 
! need to get pixel id and corresponding row id of f_p
        n = (srt_pxl-1)*vxl_cnt 
!        nf = stp_pxl*vxl_cnt
        nf = n
        do 100 k = srt_pxl, stp_pxl ! which pixels
           x_pos = MOD(k,fpa_x)
           y_pos = (k-x_pos)/fpa_x + 1 ! get y position
           if (x_pos.eq.0) then ! final check on x position
              x_pos = fpa_x
           endif
           ! is interpreted as row sparse column successive
           call spsMAT_shift(hbasis,x_pos,y_pos,fpa_x,hshift) 
           do 120 m = 1, vxl_cnt ! decompse to voxels
              call strip(hrow,hshift,m) ! gets shifted row/col
              n = n + 1 ! increment n to correct element
              tmpsum = 0.
              do 110 l = hrow(m,1),hrow(m,2) ! indices of h_mn
              if (g_k(INT(hshift(l,2))).gt.0) then
                 tmpsum = tmpsum + 
     ~        hshift(l,1)*g_o(INT(hshift(l,2)))/g_k(INT(hshift(l,2)))
              else !take it all
                 tmpsum = tmpsum +
     ~           hshift(l,1)*g_o(INT(hshift(l,2)))
              endif
110           continue
              f_k(n)=f_k(n)*(tmpsum/h_mod(m))
120        continue 
100     continue

!======================================================================
!   Send basic construct back to master task
!======================================================================

       mtype = FROM_WORKER
       call MPI_SEND(f_k(nf+1:n), (n-nf), MPI_DOUBLE_PRECISION,
     .     MASTER, mtype, MPI_COMM_WORLD, ierr)
       endif
       call MPI_FINALIZE(ierr)

       end subroutine 
