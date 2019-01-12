!======================================================================
! FILENAME: sps_BP.f
!
! DESCRIPTION:
!    This subroutine gives the initial f_k back projection utilizing
!    the generalized shift invariant inverse routine with the observed
!    data g_obs.dat
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
       subroutine sps_BP(
     &            f_k, ! returns an estimate of f_o
     &            g_k,
     &            hpinv, ! basis matrix of H pseudo inverse
     &            h_m,h_n,vxl_cnt,fpa_x,taskid)

       use com_methods

       REAL, DIMENSION(h_n), INTENT(OUT) :: f_k
       REAL, DIMENSION(h_m), INTENT(IN) :: g_k
       REAL, DIMENSION(max_hpinv,2), INTENT(IN) :: hpinv
         
       integer, parameter :: MASTER = 0
       integer, parameter :: FROM_MASTER = 1
       integer, parameter :: FROM_WORKER = 2

       integer status(MPI_STATUS_SIZE), ierr,i,j,k,l,m,n
       integer numtasks,taskid,numworkers,source,dest,nbytes,mtype
       integer srt_pxl, stp_pxl, extra, avg_pxl,yi,ye
       integer x_pos, y_pos
       real, dimension(h_n) :: f_p
       real, dimension(max_hpinv,2) :: hpinv_s
       integer, DIMENSION(vxl_cnt,2) :: hrow
       real tmpsum
!======================================================================
! f_k is shared amongst all, g_obs is broadcast to all
! preallocated
!======================================================================

!*************************** master task *****************************
       if (taskid .eq. MASTER) then
       f_k = 0.
       f_p = 0.
       extra = mod(numworkers,h_n/vxl_cnt)
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
       call MPI_SEND(g_k, h_m, MPI_DOUBLE_PRECISION,
     .    dest, mtype, MPI_COMM_WORLD, ierr)
       call MPI_SEND(hpinv, SIZE(hpinv), MPI_DOUBLE_PRECISION,
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
!     Receive f_k results from worker tasks
!     f_k sections
!======================================================================
       mtype = FROM_WORKER
       do 60 i=1, numworkers
       source = i
       call MPI_RECV(f_p, h_n, MPI_DOUBLE_PRECISION,
     .    source, mtype, MPI_COMM_WORLD, status, ierr)
       call MPI_RECV(srt_pxl, 1, MPI_INTEGER, source, mtype,
     .    MPI_COMM_WORLD, status, ierr)
       call MPI_RECV(stp_pxl, 1, MPI_INTEGER, source, mtype,
     .    MPI_COMM_WORLD, status, ierr)
!======================================================================
!gather elements of f_k   
!======================================================================
        yi =(srt_pxl-1)*vxl_cnt +1
        ye =stp_pxl*vxl_cnt
         f_k(yi:ye)=f_p(yi:ye)
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
       call MPI_RECV(g_k, h_m, MPI_DOUBLE_PRECISION, 
     .     MASTER, mtype, MPI_COMM_WORLD, status, ierr)
       call MPI_RECV(hpinv, SIZE(hpinv), MPI_DOUBLE_PRECISION,
     .     MASTER, mtype, MPI_COMM_WORLD, status, ierr)
       call MPI_RECV(srt_pxl, 1, MPI_INTEGER, MASTER,
     .     mtype, MPI_COMM_WORLD, status, ierr)
       call MPI_RECV(stp_pxl, 1, MPI_INTEGER, 
     .     MASTER, mtype, MPI_COMM_WORLD, status, ierr)
       call MPI_RECV(fpa_x, 1, MPI_INTEGER,
     .     MASTER, mtype, MPI_COMM_WORLD, status, ierr)

!======================================================================
!
!      F = H.G           Do matrix multiply 
!
!======================================================================
! 
! need to get pixel id and corresponding row id of f_p
        f_p = 0.
        n = (srt_pxl-1)*vxl_cnt
        do 100 k = srt_pxl, stp_pxl ! which pixels
           x_pos = MOD (k,fpa_x)
           y_pos = (k-x_pos)/fpa_x + 1 ! get y position
           if (x_pos.eq.0) then ! final check on x position
           x_pos = fpa_x
           endif
           ! should still be col sparse successive row
           call spsMAT_shift(hpinv,x_pos,y_pos,fpa_x,hpinv_s) 
           do 120 m = 1, vxl_cnt ! decompse to voxels
              call strip(hrow,hpinv_s,m) ! gets shifted row/col
              n = n + 1
              tmpsum = 0.0
              do 110 l = hrow(m,1),hrow(m,2) ! indices of hpinv
              tmpsum = tmpsum + hpinv_s(l,1)*g_k(INT(hpinv_s(l,2)))
110           continue
           f_p(n) = tmpsum
120        continue 
100     continue

!======================================================================
!   Send basic construct back to master task
!======================================================================

       mtype = FROM_WORKER
       call MPI_SEND(f_p, h_n, MPI_DOUBLE_PRECISION,
     .     MASTER, mtype, MPI_COMM_WORLD, ierr)
       call MPI_SEND(srt_pxl, 1, MPI_INTEGER, MASTER,
     .     mtype, MPI_COMM_WORLD, ierr)
       call MPI_SEND(stp_pxl, 1, MPI_INTEGER, MASTER,
     .     mtype, MPI_COMM_WORLD, ierr)
       endif
       call MPI_FINALIZE(ierr)

       end subroutine sps_BP 
