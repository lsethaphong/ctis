!======================================================================
! FILENAME: sps_COM.f
!
! DESCRIPTION:
! ============
!
!
!
! History of Code Changes:
! ========================
!
!   February 7, 2005: Original Code
!
!   Author:           Latsavongsakda Sethaphong
!   Institution:      Molecular Physiology & Biophysics
!                     Piston Lab
!                     Vanderbilt University
!                     Nashville, TN
!   E-mail:           latsavongsakda.sethapong@vanderbilt.edu
!----------------------------------------------------------------------

       module com_methods
       IMPLICIT NONE

!       integer, parameter :: r8 = SELECTED_REAL_KIND(15,307)
!       integer, parameter :: r4 = SELECTED_REAL_KIND(6,37)
!       integer, parameter :: i1 = SELECTED_INT_KIND(2)
!       integer, parameter :: i2 = SELECTED_INT_KIND(4)
!       integer, parameter :: i4 = SELECTED_INT_KIND(8)              
!=================================
! terms of design
!---------------------------------
       integer, parameter :: vxl_cnt = 41
       integer, parameter :: src_x = 8 ! orig 433
       integer, parameter :: src_y = 8 ! orig 433
       integer, parameter :: fpa_x = 4096
       integer, parameter :: fpa_y = 2048
       integer, parameter :: h_m = 4096*2048
       integer, parameter :: h_n = 8*8*41 ! orig 433*433*41
       integer, parameter :: max_hpinv = 200000
       integer, parameter :: max_hbasis = 41*127+1

       INTEGER, PARAMETER :: MASTER = 0
       INTEGER, PARAMETER :: FROM_MASTER = 1
       INTEGER, PARAMETER :: FROM_WORKER = 2
!==================================
! dispersion routine tags -- SEND
!----------------------------------
       INTEGER, PARAMETER :: FWRD_GK_S = 3
       INTEGER, PARAMETER :: BACK_FK_S = 4
       INTEGER, PARAMETER :: MAXF_FK_S = 5
       INTEGER, PARAMETER :: THSH_GK_S = 6
       INTEGER, PARAMETER :: MAXF_GK_S = 7
!=====================================
! master return routine tags -- RECEIVE
!-------------------------------------
       INTEGER, PARAMETER :: FWRD_GK_R = 8
       INTEGER, PARAMETER :: BACK_FK_R = 9
       INTEGER, PARAMETER :: MAXF_FK_R = 10
       INTEGER, PARAMETER :: THSH_GK_R = 11
       INTEGER, PARAMETER :: MAXF_GK_R = 12
!=====================================
! inversion routine tags -- SEND
!-------------------------------------
       INTEGER, PARAMETER :: PCI_AK_S = 13
       INTEGER, PARAMETER :: PCI_LK_S = 14
       INTEGER, PARAMETER :: PCI_LR_S = 15
       INTEGER, PARAMETER :: PCI_LI_S = 16
       INTEGER, PARAMETER :: PCI_MI_S = 17
!=====================================
! inversion routine tags -- RECIEVE
!-------------------------------------
       INTEGER, PARAMETER :: PCI_AK_R = 18
       INTEGER, PARAMETER :: PCI_LK_R = 19
       INTEGER, PARAMETER :: PCI_LR_R = 20
       INTEGER, PARAMETER :: PCI_LI_R = 21
       INTEGER, PARAMETER :: PCI_MI_R = 22
!=====================================
! terms of selection
!-------------------------------------
       INTEGER, PARAMETER :: limit = 20
       REAL, PARAMETER :: tol = 3.9215696e-3

       contains

!======================================================================
! It is more efficient to return the start and stop indices stored by
! voxel id than to generate a series of transposed vectors 
!
!----------------------------------------------------------------------
          subroutine strip(o_idx,imat, vxl_cnt)
!
!   GENERAL DESCRIPTION:
!   ====================
!
!
!
!   HISTORY of Code Changes:
!   ========================
!       February 13, 2005:   Original Code
!       
!       Author:              Latsavongsakda Sethaphong
!       Institution:         Molecular Physiology & Biophysics
!                            Piston Lab
!                            Vanderbilt University
!                            Nashville, TN
!       E-mail:              latsavongsakda.sethaphong@vanderbilt.edu
!----------------------------------------------------------------------
             INTEGER, INTENT(IN) ::  vxl_cnt
             INTEGER, DIMENSION(vxl_cnt+1), INTENT(OUT) :: o_idx
             REAL,DIMENSION(:,:), INTENT(IN) :: imat
             INTEGER cnt,i,j
             cnt = 1
             j = 1
             o_idx(1) = 1
             do 210 i = 2, SIZE(imat,1) ! go to the end
                if (imat(i,1).eq.0.and.imat(i,2).eq.0) then
                   goto 215
                endif
                if (imat(i,2).lt.0) then ! looking for the delimiter
                   j=j+1
                   o_idx(j) = i! stop
                endif
210          continue
215          return
          end subroutine strip
!======================================================================
!  END STRIP ROUTINE
!======================================================================
!======================================================================
! Because the inverse matrix will move in tandem with forward matrix 
! only one function is required.  The inverse is column sparse whereas
! the forward matrix is row sparse
! the first row is reserved for dimension info e.g. rows, col
!----------------------------------------------------------------------
         subroutine spsMAT_shift(imat, rshift, 
     ~ dshift, fpa_x, smat)
!   GENERAL DESCRIPTION:
!   ====================
!
!
!
!   HISTORY of Code Changes:
!   ========================
!       February 13, 2005:   Original Code
!
!       Author:              Latsavongsakda Sethaphong
!       Institution:         Molecular Physiology & Biophysics
!                            Piston Lab
!                            Vanderbilt University
!                            Nashville, TN
!       E-mail:              latsavongsakda.sethaphong@vanderbilt.edu
!----------------------------------------------------------------------

             REAL, INTENT(IN), DIMENSION(:,:) :: imat
             REAL, INTENT(OUT), DIMENSION(:,:) :: smat
             INTEGER, INTENT(IN) :: rshift
             INTEGER, INTENT(IN) :: dshift
             INTEGER, INTENT(IN) :: fpa_x
             INTEGER mdshift,i
!            ----------------------------------
             mdshift = dshift*fpa_x + rshift
             smat = imat !try this first
             do 310 i = 2, SIZE(imat,1)
                if(imat(i,1).eq.0.and.imat(i,2).eq.0) then
                   goto 315
                endif
!                if (imat(i,2).lt.0) then
!                   smat(i,2) = -(abs(imat(i,2))+mdshift)
!                else
              smat(i,2) = abs(imat(i,2)) + mdshift
!                endif
310          continue
315          return             
         end subroutine  spsMAT_shift

!=====================================================================
! END SPSMAT_SHIFT
!=====================================================================

!---------------------------------------------------------------------
         subroutine MOD_H(imat, mat_mod)
!
!   GENERAL DESCRIPTION:
!   ====================
!
!
!
!   HISTORY of Code Changes:
!   ========================
!       February 13, 2005:   Original Code
!
!       Author:              Latsavongsakda Sethaphong
!       Institution:         Molecular Physiology & Biophysics
!                            Piston Lab
!                            Vanderbilt University
!                            Nashville, TN
!       E-mail:              latsavongsakda.sethaphong@vanderbilt.edu
!----------------------------------------------------------------------
         REAL, DIMENSION(:,:), INTENT(IN) :: imat
         REAL, DIMENSION(vxl_cnt), INTENT(OUT) :: mat_mod
         INTEGER i,j
         REAL tmpsum ! may need to adjust buffer size
         tmpsum = 0.
         j = 0
         do 400 i = 2,SIZE(imat,1)
            if (imat(i,1).eq.0.and.imat(i,2).eq.0) then
               goto 405
            endif
            tmpsum=tmpsum+imat(i,1)
            if (imat(i,2).lt.0) then
               j = j+1
               mat_mod(j) = tmpsum
               tmpsum = 0. 
            endif 
400      continue
405      return
         end subroutine mod_h
!=====================================================================
!  END MOD_H
!=====================================================================

!=====================================================================
! ########################## END OF FILE #############################
!---------------------------------------------------------------------
       end module com_methods
