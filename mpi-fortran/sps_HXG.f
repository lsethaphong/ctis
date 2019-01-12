!======================================================================
! This program calculates the estimate of f after obtaining an initial
! guess thru back projection of the observed focal plane image pattern.
! 
!
!======================================================================

       program sps_HXG(hname,col_start, col_stop, g_obs, g_est,f_est)
       character*25 hname, g_obs, g_est, f_est
       integer col_start, col_stop
!
! Load up the basis matrix H
!
       nz_h_vxl = 
       do 5 j = col_start, col_stop ! n of Hmn
          do 10 i = 1,nz_h_vxl ! m of Hmn -- do all m from  basis h
             
10        end do
5      end do
       end
