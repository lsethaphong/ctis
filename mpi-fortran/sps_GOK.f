!======================================================================
! FILENAME: sps_GOK.f
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
       subroutine sps_GOK(g_o,g_k,sps_stop,taskid) 
       use com_methods
       REAL, DIMENSION(:), INTENT(IN) :: g_o
       REAL, DIMENSION(:), INTENT(IN) :: g_k
       LOGICAL, INTENT(OUT) :: sps_stop
       integer ierr, i,j,k,l,m,n
       integer numtasks,taskid,numworkers,source,dest,nbytes,mtype
       real my_diff, final_diff
       integer length, slice, i_start, i_stop, tag, item
       integer j_start,j_stop
!======================================================================
! f_k is shared amongst all, g_o is broadcast to all
!======================================================================

     
       numworkers = numtasks-1
       length = size(g_k,1)
       slice = (length/numtasks)!master also performs the function
       i_start = taskid*slice + 1
       i_stop  = (taskid+1)*slice 
!       my_diff = euclidean_diff(g_k,g_o,i_start,i_stop)

      if (taskid.eq.MASTER) then
! sending whole vectors
           do source = 1, numworkers
                   j_start = taskid*slice + 1
                   i_stop = (source+1)*slice  
                   
!          write(*,*) '   sending',(j_stop-istart),' cols to task',source
              call MPI_SEND(g_k(j_start:j_stop),slice,MPI_REAL, source, 
     .             0, MPI_COMM_WORLD, ierr)
                   
              call MPI_SEND(g_o(j_start:j_stop),slice,MPI_REAL, source, 
     .             1, MPI_COMM_WORLD, ierr)
           enddo
           final_diff = euclidean_diff(g_k(i_start:i_stop),
     .             g_o(i_start:i_stop),i_start,i_stop)
           do source = 1, numworkers
!           write(*,*) 'calling ',source
               call MPI_RECV(my_diff, 1, MPI_REAL, source, FROM_WORKER,
     .              MPI_COMM_WORLD, status, ierr)
              final_diff  = final_diff + my_diff
!           write(*,*) '   got ',my_diff
          enddo
          if (final_diff < tol) then
             sps_stop = .true.
          else
             sps_stop = .false.
          endif
      else
! recieve g_k and g_o from the master 
          call MPI_RECV(g_k(i_start:i_stop),slice, MPI_REAL, 0, 
     .         0, MPI_COMM_WORLD, status, ierr)
          call MPI_RECV(g_o(i_start:i_stop),slice, MPI_REAL, 0,
     .         1, MPI_COMM_WORLD, status, ierr)
!          dest = 0 ! all go to Master
! returning only the partial sum
          my_diff = euclidean_diff(g_k,g_o,i_start,i_stop)
          call MPI_SEND(my_diff, 1, MPI_REAL, 0,
     .         FROM_WORKER, MPI_COMM_WORLD, ierr)
      endif
      end subroutine sps_GOK


       real function euclidean_diff(g_i,g_s,i,s)
       IMPLICIT NONE
       REAL, DIMENSION(:), INTENT(IN) :: g_i, g_s
       REAL sqr_dif, tmp_dif
       INTEGER i,s,k ! i start, j stop
       sqr_dif =0.
       tmp_dif = 0.
       do 400 k=i,s
       tmp_dif = g_i(k)-g_s(k)
       sqr_dif = sqr_dif + tmp_dif*tmp_dif
400    continue
       euclidean_diff = sqr_dif 
       end function euclidean_diff
