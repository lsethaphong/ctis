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
       program sps_GOK 

       REAL, DIMENSION(8) :: g_o
       REAL, DIMENSION(8) :: g_k
       LOGICAL :: sps_stop

       integer, parameter :: MASTER = 0
       integer, parameter :: FROM_MASTER = 1
       integer, parameter :: FROM_WORKER = 2
       real, parameter :: tol = 1.e-5

       include 'mpif.h'
       integer status(MPI_STATUS_SIZE), ierr, i,j,k,l,m,n
       integer numtasks,taskid,numworkers,source,dest,nbytes,mtype
       real my_diff, final_diff
       integer length, slice, i_start, i_stop, tag, item
       integer j_start, j_stop
!======================================================================
! f_k is shared amongst all, g_o is broadcast to all
!======================================================================

       call MPI_INIT(ierr) ! initialize mpi system
       call MPI_COMM_RANK(MPI_COMM_WORLD, taskid, ierr)
       call MPI_COMM_SIZE(MPI_COMM_WORLD, numtasks, ierr)
     
       numworkers = numtasks-1
       length = size(g_k,1)
       slice = (length/numtasks)!master also performs the function
       i_start = taskid*slice + 1
       i_stop  = (taskid+1)*slice 
!       my_diff = euclidean_diff(g_k,g_o,i_start,i_stop)

      if (taskid.eq.MASTER) then
       
       do i=1,8
          g_k(i) = i
          g_o(i) = i
       enddo
!       write(*,*) ' written g_k & g_o'
!       write(*,*) ' g_k ', g_k
!       write(*,*) ' g_o ', g_o
! sending whole vectors 
           if (numworkers > 0) then
           do source = 1, numworkers
! partition the data
              j_start = source*slice + 1
              j_stop = (source +1)*slice
                   item = 0
           write(*,*) '   sending',slice,' to task', source
              call MPI_SEND(g_k(j_start:j_stop),slice,MPI_REAL, source, 
     .             item, MPI_COMM_WORLD, ierr)
                   item = 1
              call MPI_SEND(g_o(j_start:j_stop),slice,MPI_REAL, source, 
     .             item, MPI_COMM_WORLD, ierr)
           enddo
           endif
           final_diff = euclidean_diff(g_k(i_start:i_stop),
     .             g_o(i_start:i_stop),i_start,i_stop)
           write(*,*) taskid, ' master ', final_diff
           if (numworkers > 0) then
           do source = 1, numworkers
           write(*,*) 'calling ',source
               call MPI_RECV(my_diff, 1, MPI_REAL, source, FROM_WORKER,
     .              MPI_COMM_WORLD, status, ierr)
              final_diff  = final_diff + my_diff
           write(*,*) '   got ',my_diff
          enddo
          endif
          if (final_diff < tol) then
             sps_stop = .true.
          else
             sps_stop = .false.
          endif
          write(*,*) ' final answer', sps_stop
      else
! recieve g_k and g_o from the master 
          call MPI_RECV(g_k(i_start:i_stop),slice, MPI_REAL, 0, 
     .         0, MPI_COMM_WORLD, status, ierr)
          write(*,*) 'recieved', g_k(i_start:i_stop)
          call MPI_RECV(g_o(i_start:i_stop),slice, MPI_REAL, 0,
     .         1, MPI_COMM_WORLD, status, ierr)
          write(*,*) 'recieved', g_o(i_start:i_stop)
!          dest = 0 ! all go to Master
! returning only the partial sum
          write(*,*) taskid, ' ', slice, 'is slice'
          write(*,*) taskid, ' got gk: ',g_k(i_start:i_stop)
          my_diff = euclidean_diff(g_k,g_o,i_start,i_stop)
          write(*,*) taskid, ' worker ', my_diff
          call MPI_SEND(my_diff, 1, MPI_REAL, 0,
     .         FROM_WORKER, MPI_COMM_WORLD, ierr)
      endif
          call MPI_FINALIZE(ierr)
      end program


       real function euclidean_diff(g_i,g_s,i,s)
       IMPLICIT NONE
       REAL, DIMENSION(8), INTENT(IN) :: g_i, g_s
       integer, intent(in) :: i
       integer, intent(in) :: s
       REAL sqr_dif, tmp_dif
       integer k ! i start, j stop
       sqr_dif =0.
       tmp_dif = 0.
       write(*,*) ' working here, ',i,',',s
       write(*,*) ' is ', g_i(i:s)
       do 400 k=i,s
       tmp_dif = g_i(k)-g_s(k)
       sqr_dif = sqr_dif + tmp_dif*tmp_dif
400    continue
       euclidean_diff = sqr_dif 
       return
       end function euclidean_diff


!       real function square(x)
!       implicit none
!       real x
!       square = x*x
!       end function square
