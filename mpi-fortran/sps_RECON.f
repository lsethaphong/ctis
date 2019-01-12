!======================================================================
! FILENAME:
!    sps_RECON.f
! DESCRIPTION:
!    This code calls up the parallelized reconstruction subroutines 
!    thereby implementing multiplicative algebraic reconstrunction
!    on a file desiganted g_obs.dat, the experimental data
!    to produce an object plane source  
!
! HISTORY:
!   February 7, 2005: Original Code
!   March 8, 2005   : Expanded routines to include SPS_MAXF, SPS_FP
!   
!   Author:           Latsavongsakda Sethaphong
!   Institution:      Molecular Physiology & Biophysics
!                     Piston Lab
!                     Vanderbilt University
!                     Nashville, TN
!                     latsavongsakda.sethapong@vanderbilt.edu
!----------------------------------------------------------------------
       include 'sps_COM.f'
       program SPS_RECON
       use com_methods
       IMPLICIT NONE
       REAL, DIMENSION(h_m) :: g_k
       REAL, DIMENSION(h_m) :: g_o
       REAL, DIMENSION(h_n) :: f_k
       include 'mpif.h'
! Basis matrix projection and inverse
       integer status(MPI_STATUS_SIZE), ierr, gonow
       REAL, DIMENSION(max_hbasis,2) :: hbasis
       REAL, DIMENSION(max_hpinv,2) :: hpinv
       REAL, DIMENSION(vxl_cnt) :: h_mod ! column algebraic modulus
       INTEGER, DIMENSION(vxl_cnt+1) :: hrow
       INTEGER i,j,k,taskid,numtasks,numworkers
       real e,etime,t(2),maxfkn,mu_err      
!======================================
! INITIALIZE MPI AND GET RANK, UNIVERSE
!--------------------------------------
       CALL MPI_INIT(ierr)
       CALL MPI_COMM_RANK(MPI_COMM_WORLD, taskid, ierr)
       CALL MPI_COMM_SIZE(MPI_COMM_WORLD, numtasks, ierr)
!=======================================
! UPLOAD BASIS MATRIX AND PSEUDO INVERSE
! IN COLUMN MAJOR ORDER
!---------------------------------------
!       write(*,*) taskid, 'ready'

       open(unit=12, file='hbasis.dat',status='old')
       i = 1
       do while(.TRUE.)
          read(12,end=35,fmt=*) hbasis(i,1), hbasis(i,2)
          i = i+1
       enddo
35     close(12)
!=========================
! UPLOAD THE INVERSE BASIS
! COLUMN MAJOR ORDER
!-------------------------
       open(unit=13, file='hpinv.dat',status='old')
       i = 1
       do while(.TRUE.)
          read(13,end=45,fmt=*) hpinv(i,1), hpinv(i,2)
          i=i+1
       enddo
45     close(13)
 
       numworkers = numtasks - 1

!=====================================================================
! BEGIN MAIN TASKID SPECIFIC FUNCTIONS ###############################
!---------------------------------------------------------------------

       if (taskid.eq.MASTER) then
          e=etime(t)
       endif

       open(unit=11, file='g_o.dat', status='old')
       i = 1
       do 50 while(.TRUE.) ! length of data
           read(11,end=55,fmt=*) g_o(i)
           i = i+1
50     enddo
55     close(11)

!--------------------------------------------
! Synchronize Pause
!--------------------------------------------
!       call MPI_BARRIER(MPI_COMM_WORLD, ierr)
!--------------------------------------------

!============================================
! ONLY THE MASTER READS/DISTRIBUTES THE DATA
!--------------------------------------------
! Check dimensions for the rest of the algorithm
! GET F_K from G_O       
60     call SPS_BP(f_k,g_o,hpinv,
     .      taskid,numworkers)
!============================================
! get modulus of columns in hbasis and store 
! in h_mod from com_methods
!--------------------------------------------
       call strip(hrow,hbasis,vxl_cnt)
       call MOD_H(hbasis,h_mod)
!--------------------------------------------
! THE MART ALGORITHM
!============================================
       do i = 1, limit
!--------------------------------------------
! Forward Projection to get g_k
!       if (taskid.eq.MASTER) then
!          write(*,*) 'MAXF iteration ',i
!       endif
!--------------------------------------------
       call SPS_FP(f_k,g_k,hbasis,hrow,
     .      taskid,numworkers) 

!--------------------------------------------
! Calculate euclidian difference
! r(e) = sigma(Gm - Mm)**2)
!--------------------------------------------
!--------------------------------------------
! Implement positive saturation thresholding
!--------------------------------------------
       call SPS_SAT(g_k,g_o,mu_err,taskid,numworkers)
!--------------------------------------------
! Correct f_k by scaling of h_mod and shgom
       call SPS_MAXF(f_k,h_mod,hbasis,g_k,hrow,
     .      taskid,numworkers)
!--------------------------------------------

!--------------------------------------------
! Synchronize Pause
!--------------------------------------------
       call MPI_BARRIER(MPI_COMM_WORLD, ierr)
!--------------------------------------------
! END OF MAXF LOOP
!============================================       
65     enddo

!============================================
! full data is kept with the MASTER
!--------------------------------------------
       if (taskid.eq.MASTER) then
!         maxfkn = MAXVAL(f_k)
         open(unit=16, file='f_f.dat', status='unknown')
         do 70 i=1,h_n
!               write(16,'(F5.3)') f_k(i)/maxfkn
               write(16,'(F5.3)') f_k(i)
70       continue
         close(16)
         maxfkn = MAXVAL(g_k)
         open(unit=17, file='g_k.dat',status='unknown')
         do 80 i =1,h_m
            write(17,'(F5.3)') g_k(i)/maxfkn
80       continue
         close(17)
        e = etime(t)
       write(*,*) '---------------------------------------------------'
       print *, 'elapsed:', e, 'user:', t(1), 'sys:', t(2)
       write(*,*) 'COMPLETED PROGRAM'
       endif

!---------------------------------------------
!=====================================================================
!       END OF MAIN PROGRAM
!---------------------------------------------------------------------
       call MPI_FINALIZE(ierr)
100    end

!======================================================================
! ROUTINE NAME: sps_BP
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
       subroutine SPS_BP(f_k,g_k,hpinv,taskid,numworkers)
       use com_methods
       IMPLICIT NONE
       include 'mpif.h'
       real, dimension(h_n), intent(OUT) :: f_k
       real, dimension(h_m), intent(IN) :: g_k
       real, dimension(max_hpinv,2), intent(IN) :: hpinv
       integer, intent(IN) :: taskid
       integer, intent(IN) :: numworkers
       integer status(MPI_STATUS_SIZE)
       integer i_start, i_stop, slice,s,z,i,mtype
       integer srt_pxl,stp_pxl,avg_pxl,extra,ierr
       
!determine extra pixels
       extra = mod((numworkers+1),h_n/vxl_cnt)
       avg_pxl = (h_n/vxl_cnt-extra)/(numworkers+1)
       extra = h_n/vxl_cnt-(avg_pxl*(numworkers+1))

!taskid starts at zero
       if (taskid.eq.MASTER) then
          write(*,*) 'Backward Projecting'
          srt_pxl = 1
          stp_pxl = (taskid+1)*avg_pxl + extra
       else
          srt_pxl = taskid*avg_pxl + extra + 1 
          stp_pxl = (taskid+1)*avg_pxl + extra
       endif

       f_k = 0.

       call MPI_BCAST(g_k,h_m, MPI_REAL,
     .    MASTER, MPI_COMM_WORLD, ierr)
 
       call CORE_BP(f_k,g_k,srt_pxl,stp_pxl,hpinv)

       if (taskid.eq.MASTER) then
! call back based on vxl_cnt elements 
          mtype = BACK_FK_R      
          do i =1,numworkers
             s = (i*avg_pxl + extra)*vxl_cnt + 1 ! start
             z = ((i+1)*avg_pxl + extra)*vxl_cnt ! stop
! storing in local buffer f_k
          call MPI_RECV(f_k(s:z),z-s+1, MPI_REAL,
     .     i, mtype, MPI_COMM_WORLD, status, ierr) 
!          write(*,*) 'sent ',taskid,s,z
          enddo
       else
!======================================================================
!   Send basic construct back to master task
!======================================================================
       s = (srt_pxl-1)*vxl_cnt + 1
       z = stp_pxl*vxl_cnt
       mtype = BACK_FK_R
       call MPI_SEND(f_k(s:z),z-s+1, MPI_REAL,
     .     MASTER, mtype, MPI_COMM_WORLD, ierr)
!       write(*,*) taskid, 'sent f_k',s,z
       endif
       end subroutine SPS_BP

       subroutine CORE_BP(f_k,g_k,srt_pxl,stp_pxl,hpinv)
       use com_methods
! need to get pixel id and corresponding row id of f_p
       IMPLICIT NONE
       real, dimension(h_n), intent(OUT) :: f_k
       real, dimension(h_m), intent(IN) :: g_k
       integer, intent(IN) :: srt_pxl
       integer, intent(IN) :: stp_pxl
       real, dimension(max_hpinv,2), intent(IN) :: hpinv
!       real, dimension(max_hpinv,2) :: hpinv_s
       integer, dimension(vxl_cnt+1) :: hrow
       real tmpsum     
       integer k,l,m,n,x_pos,y_pos,mdshift,idx

! use src_x for determing pixels
       call strip(hrow,hpinv,vxl_cnt) ! gets shifted row/col
        do 200 k = srt_pxl, stp_pxl ! which pixels
           x_pos = MOD(k,src_x) ! src_x is defined in com_methods
           if (x_pos.eq.0) then ! final check on x position
              x_pos = src_x
           endif
           y_pos = (k-x_pos)/src_x
           n = (k-1)*vxl_cnt ! gets the starting voxel
!           call spsMAT_shift(hpinv,x_pos-1,y_pos-1,fpa_x,hpinv_s)
           mdshift = y_pos*fpa_x+x_pos-1
           do 220 m = 1, vxl_cnt ! decompse to voxels
             tmpsum=0.
              do 210 l = (hrow(m)+1),hrow(m+1) ! indices of hpinv
                 idx=INT(ABS(hpinv(l,2)))+mdshift
                 tmpsum=tmpsum+hpinv(l,1)*g_k(idx)
210           continue
              if (tmpsum.gt.1e-6) then ! omit negative values
                 f_k(n+m)=tmpsum
              endif
!              if (tmpsum.ge.1.000) then
!                 write(*,*) 'inverse confirm,',tmpsum,'at',n+m
!              endif
220        continue
200     continue

225     endsubroutine CORE_BP
!---------------------------------------------------------------------
! END OF BACK PROJECTION ROUTINES
!=====================================================================

!=====================================================================
! ROUTINE NAME: sps_FP
!
! DESCRIPTION:
!    This subroutine gives the initial g_k forward projection utilizing
!    the generalized shift invariant basis matrix H
!
! HISTORY:
!   February 7, 2005: Original Code
!   March, 3, 2005  : Adapted for MPI in large format
!
!   Author:           Latsavongsakda Sethaphong
!   Institution:      Molecular Physiology & Biophysics
!                     Piston Lab
!                     Vanderbilt University
!                     Nashville, TN
!                     latsavongsakda.sethapong@vanderbilt.edu
!---------------------------------------------------------------------
       subroutine SPS_FP(f_k,g_k,hbasis,hrow,taskid,numworkers)
       use com_methods
       IMPLICIT NONE
       include 'mpif.h'
       real, dimension(h_n), intent(IN) :: f_k
       real, dimension(h_m), intent(OUT) :: g_k
       real, dimension(max_hbasis,2), intent(IN) :: hbasis
       integer, dimension(vxl_cnt+1) :: hrow
       integer, intent(IN) :: taskid
       integer, intent(IN) :: numworkers
       real, dimension(h_m) :: g_p ! temporary register
       integer status(MPI_STATUS_SIZE)
       integer s,z,i,mtype,j
       integer srt_pxl,stp_pxl,avg_pxl,extra,ierr

! determine extra pixels
       extra = mod((numworkers+1),h_n/vxl_cnt)
       avg_pxl = (h_n/vxl_cnt-extra)/(numworkers+1)
       extra = h_n/vxl_cnt-(avg_pxl*(numworkers+1))
! taskid starts at zero
       if (taskid.eq.MASTER) then
          write(*,*) 'Forward Projection'
          srt_pxl = 1
          stp_pxl = (taskid+1)*avg_pxl + extra
       else
          srt_pxl = taskid*avg_pxl + 1 + extra
          stp_pxl = (taskid+1)*avg_pxl + extra
       endif
! clearing the buffer
       g_k =0.
       g_p =0.
       call MPI_BCAST(f_k,h_n,MPI_REAL,MASTER,MPI_COMM_WORLD,ierr)
       call CORE_FP(f_k,g_k,srt_pxl,stp_pxl,hbasis,hrow)
!======================================================================
!   Send basic construct back to master task
!======================================================================
       mtype = FWRD_GK_R
       if (taskid.eq.MASTER) then
          do i = 1,numworkers
             call MPI_RECV(g_p,h_m,MPI_REAL,i,
     .            mtype,MPI_COMM_WORLD,status,ierr)
             do j=1,h_m  ! saturation mode
                g_k(j)=MIN(g_k(j)+g_p(j),1.000)
!                g_k(j)=g_k(j)+g_p(j)
             enddo
          enddo
       else
          call MPI_SEND(g_k,h_m,MPI_REAL,
     .         MASTER,mtype,MPI_COMM_WORLD,ierr)
       endif
       call MPI_BARRIER(MPI_COMM_WORLD,ierr)
       end subroutine SPS_FP

       subroutine CORE_FP(f_k,g_k,srt_pxl,stp_pxl,hbasis,hrow)
       use com_methods
! need to get pixel id and corresponding row id of f_p
       IMPLICIT NONE
       real, dimension(h_n), intent(IN) :: f_k
       real, dimension(h_m), intent(OUT) :: g_k
       integer, intent(IN) :: srt_pxl
       integer, intent(IN) :: stp_pxl
       real, dimension(max_hbasis,2), intent(IN) :: hbasis
       real, dimension(max_hbasis,2) :: hshift
       integer, dimension(vxl_cnt+1), intent(IN) :: hrow
       integer k,l,m,n,x_pos,y_pos,idx,mdshift
! switch to source dimensions
! need to get pixel id and corresponding row id of f_p
        do 300 k = srt_pxl, stp_pxl ! which pixels
           x_pos = MOD(k,src_x)
           if (x_pos.eq.0) then ! final check on x position
              x_pos = src_x
           endif
           y_pos = (k-x_pos)/src_x! get y position
           mdshift = y_pos*fpa_x+x_pos-1
           n = (k-1)*vxl_cnt ! calculate last index of f_k
           do 320 m = 1, vxl_cnt ! decompse to voxels
!              if (f_k(n+m).gt.tol) then ! skip pts of zero probability
                 do 310 l = (hrow(m)+1),hrow(m+1) ! indices of hbasis per voxel
                 idx = INT(ABS(hbasis(l,2)))+mdshift 
                 g_k(idx) = g_k(idx)+f_k(n+m)*hbasis(l,1)
310              continue
!              endif
320        continue
300     continue
       end subroutine CORE_FP

!---------------------------------------------------------------------
! END OF FORWARD PROJECTION ROUTINES
!=====================================================================

!=====================================================================
! BEGIN MAXIMIZATION CALCULATIONS
!---------------------------------------------------------------------

!======================================================================
! ROUTINE NAME: SPS_MAXF
!
! DESCRIPTION: Gives the corrected f_k given g_o and the forward
!    projection g_k from f_k; the corrected f_k is returened.
!
!
! HISTORY:
!   February 7, 2005: Original Code
!   March 7, 2005   : Broken into a core and calling routine   
!
!   Author:           Latsavongsakda Sethaphong
!   Institution:      Molecular Physiology & Biophysics
!                     Piston Lab
!                     Vanderbilt University
!                     Nashville, TN
!                     latsavongsakda.sethapong@vanderbilt.edu
!----------------------------------------------------------------------

       subroutine SPS_MAXF(f_k,h_mod,hbasis,g_w,hrow,taskid,numworkers)
       use com_methods
       IMPLICIT NONE
       include 'mpif.h'
       REAL, DIMENSION(h_n), INTENT(INOUT)  :: f_k
       REAL, DIMENSION(vxl_cnt), INTENT(IN) :: h_mod
       REAL, DIMENSION(max_hbasis,2), INTENT(IN) :: hbasis
!       REAL, DIMENSION(h_m), INTENT(IN) :: g_o
       INTEGER, DIMENSION(vxl_cnt+1), INTENT(IN) :: hrow
       REAL, DIMENSION(h_m), INTENT(IN) :: g_w ! weighted 
       INTEGER, INTENT(IN) :: taskid
       INTEGER, INTENT(IN) :: numworkers
       integer status(MPI_STATUS_SIZE)
       integer s,z,i,mtype
       integer srt_pxl,stp_pxl,avg_pxl,extra,ierr

!determine extra pixels
       extra = mod((numworkers+1),h_n/vxl_cnt)
       avg_pxl = (h_n/vxl_cnt-extra)/(numworkers+1)
       extra = h_n/vxl_cnt-(avg_pxl*(numworkers+1))

!taskid starts at zero
          if (taskid.eq.MASTER) then
             write(*,*) 'Starting MAXF'
             srt_pxl = 1
             stp_pxl = (taskid+1)*avg_pxl + extra
          else
             srt_pxl = taskid*avg_pxl + extra + 1
             stp_pxl = (taskid+1)*avg_pxl + extra
          endif

          s = (srt_pxl-1)*vxl_cnt + 1
          z = stp_pxl*vxl_cnt

       call MPI_BCAST(g_w,h_m,MPI_REAL,MASTER,MPI_COMM_WORLD,ierr)
       call CORE_MAXF(f_k,h_mod,srt_pxl,stp_pxl,hbasis,g_w,hrow)

       if (taskid.eq.MASTER) then
! call back based on vxl_cnt elements
          mtype = MAXF_FK_R
          do i =1,numworkers
             s = (i*avg_pxl + extra)*vxl_cnt + 1 ! start
             z = ((i+1)*avg_pxl + extra)*vxl_cnt ! stop
! storing in local buffer f_k
             call MPI_RECV(f_k(s:z),z-s+1,MPI_REAL,
     .            i,mtype,MPI_COMM_WORLD,status,ierr)
!             write(*,*) 'received',taskid,s,z
          enddo
!======================================================================
!   Send basic construct back to master task
!======================================================================
       else
          mtype = MAXF_FK_R
          call MPI_SEND(f_k(s:z),z-s+1,MPI_REAL,
     .         MASTER,mtype,MPI_COMM_WORLD,ierr)
       endif
       end subroutine SPS_MAXF

       subroutine CORE_MAXF(f_k,h_mod,srt_pxl,stp_pxl,hbasis,g_w,hrow)
       use com_methods
       IMPLICIT NONE
       REAL, DIMENSION(h_n), INTENT(INOUT) :: f_k
       REAL, DIMENSION(vxl_cnt), INTENT(IN) :: h_mod
       REAL, DIMENSION(h_m), INTENT(IN) :: g_w
       INTEGER, INTENT(IN) :: srt_pxl
       INTEGER, INTENT(IN) :: stp_pxl
       REAL, DIMENSION(max_hbasis,2), INTENT(IN) :: hbasis
       integer, dimension(vxl_cnt+1), INTENT(IN) :: hrow
       integer k,l,m,n,x_pos,y_pos,idx,mdshift
       real tmpsum, tmpfkn

! use source dimensions to get pixel id's
! need to get pixel id and corresponding row id of f_p
        tmpfkn=0.
!        call strip(hrow,hbasis,vxl_cnt)
        do 400 k = srt_pxl, stp_pxl ! which pixels
           x_pos = MOD(k,src_x)
           if (x_pos.eq.0) then ! final check on x position
              x_pos = src_x
           endif
           y_pos = (k-x_pos)/src_x! get y position
           n = (k-1)*vxl_cnt
           ! is interpreted as row sparse column successive
           mdshift = y_pos*fpa_x+x_pos-1
           do 420 m = 1, vxl_cnt ! decompse to voxels
              tmpsum = 0.
              if (f_k(n+m).gt.1e-6) then ! thresholding of f_k
              tmpfkn = f_k(n+m)/h_mod(m) ! normalize f_k
                 do 410 l =(hrow(m)+1),hrow(m+1)  ! indices of h_mn
                    idx = INT(ABS(hbasis(l,2)))+mdshift
                    tmpsum = tmpsum+hbasis(l,1)*g_w(idx) ! weighted
410              continue
              endif
              f_k(n+m)=tmpsum*tmpfkn
420        continue
400     continue
       end subroutine CORE_MAXF
!---------------------------------------------------------------------
! END OF MULTIPLICATIVE ALGEBRAIC RECONSTRUCTION
!=====================================================================

!=====================================================================
! BEGIN THRESHOLDING AND CORRECTION FOR G_K TO G_D
!---------------------------------------------------------------------

!======================================================================
! ROUTINE NAME: SPS_SAT
!
! DESCRIPTION: Gives the weighted g_k after thresholding and
!    scaling with g_o
!    
! HISTORY:
!   March 14, 2005   : Original Code
!
!   Author:           Latsavongsakda Sethaphong
!   Institution:      Molecular Physiology & Biophysics
!                     Piston Lab
!                     Vanderbilt University
!                     Nashville, TN
!                     latsavongsakda.sethapong@vanderbilt.edu
!----------------------------------------------------------------------

       subroutine SPS_SAT(g_k,g_o,mu_err,taskid,numworkers)
       use com_methods
       IMPLICIT NONE
       include 'mpif.h'
       REAL, DIMENSION(h_m), INTENT(INOUT):: g_k
       REAL, DIMENSION(h_m), INTENT(IN) :: g_o ! observed 
       INTEGER, INTENT(IN) :: taskid
       INTEGER, INTENT(IN) :: numworkers
       REAL, INTENT(OUT) :: mu_err
       integer status(MPI_STATUS_SIZE)
       integer srt_idx,stp_idx,avg_idx,extra,ierr
       integer i,s,z,mtype
       real mu_now
!determine extra indices
       extra = mod((numworkers+1),h_m)
       avg_idx = (h_m-extra)/(numworkers+1)
       extra = h_m-(avg_idx*(numworkers+1))

!taskid starts at zero
       if (taskid.eq.MASTER) then
          write(*,*) 'Positive Saturation Thresholding'
          srt_idx = 1
          stp_idx = (taskid+1)*avg_idx + extra
       else
          srt_idx = taskid*avg_idx + extra + 1
          stp_idx = (taskid+1)*avg_idx + extra
       endif
!---------------------------------------------------------------------
! GET THE COMMON G_K FORWARD PROJECTION
!---------------------------------------------------------------------
       s = srt_idx
       z = stp_idx

       call MPI_BCAST(g_k,h_m,MPI_REAL,MASTER,MPI_COMM_WORLD,ierr)
       call CORE_SAT(g_k(s:z),g_o(s:z),s,z,mu_now)
       call MPI_REDUCE(mu_now,mu_err,1,MPI_REAL,MPI_SUM,
     .      MASTER,MPI_COMM_WORLD,ierr)

       if (taskid.eq.MASTER) then
! call back based on vxl_cnt elements
          mtype = THSH_GK_R
          do i =1,numworkers
          s = i*avg_idx + extra + 1 ! start
          z = (i+1)*avg_idx + extra ! stop
! storing in local buffer f_k
          call MPI_RECV(g_k(s:z),z-s+1,MPI_REAL,
     .         i,mtype,MPI_COMM_WORLD,status,ierr)
       enddo
          write(*,*) 'r(e)=',mu_err
!---------------------------------------------------------------------
       else
!======================================================================
!   Send basic construct back to master task
!======================================================================
          mtype = THSH_GK_R
          call MPI_SEND(g_k(s:z),z-s+1, MPI_REAL,
     .         MASTER, mtype,MPI_COMM_WORLD,ierr)
       endif
       end subroutine SPS_SAT

       subroutine CORE_SAT(g_ws,g_os,si,st,mu_err)
       use com_methods
       INTEGER, INTENT(IN) :: si,st
       REAL, DIMENSION(si:st), INTENT(INOUT) :: g_ws
       REAL, DIMENSION(si:st), INTENT(IN) :: g_os
       INTEGER i
       REAL, INTENT(OUT) :: mu_err
       mu_err = 0.
       do i=si,st
          mu_err = mu_err + (g_os(i)-g_ws(i))**2
          if (g_ws(i).lt.1e-6) then
             g_ws(i) = g_os(i)
          else ! natural weighting
             g_ws(i)= g_os(i)/g_ws(i)
          endif
       enddo
       end subroutine CORE_SAT
!---------------------------------------------------------------------
! END OF THRESHOLDING
!=====================================================================
