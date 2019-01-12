!======================================================================
! FILENAME: sps_GOK_TB.f
!
! DESCRIPTION:
!    Test bench for sps_GOK.f subroutine
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
       include 'sps_GOK.f'
       program sps_GOK_TB
       parameter (h_n = 16)
       REAL, DIMENSION(h_n) :: g_k, g_o
       logical sps_stop
       integer i

!       sps_stop = .FALSE.
       do i=1,16
          g_k(i) = i
          g_o(i) = 2*i+1
       enddo
       write(*,*) ' written g_k & g_o'
       call sps_GOK(g_k,g_o, sps_stop)
       open(unit = 11, file='report.dat', status='new')
       write(11,*) sps_stop
       close(11)
       end
