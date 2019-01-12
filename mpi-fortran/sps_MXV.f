!======================================================================
! This function gets the object estimated for voxel elements
! corresponding to row_start and ending with row_stop integer
! arguments.  vecfile is the name of the source vector, uncompressed
! port file is the fragment file that defines the output vector
! The system matrix is row sparse, and is indexed by column.
! The forward projection is handled such that the resulting 
! g estimate is a sum of all the individual f estimat partitions
!
! The master program will sum up the composite estimates from all the 
! slave/ancillary programs to produce the final projection 
! Hence, this operation is not entirely multiplicative in nature
! For all non zero elements of the f estimate, a scaled value projection
! will be generated into a temporary summing vector. 
!
! 2/02/2005
! Latsavongsakda Sethaphong
! Piston Lab VUMC
!
!======================================================================
       subroutine sps_MXV(row_start, row_stop, hfile, 
     ~          vecfile,veclength, portfile)
       implicit integer(i-n)
       implicit real(a-h,o-z)
       character*25 hfile, vecfile, portfile
       integer  row_start, row_stop, veclength
       dimension vec(veclength)
! Read in the basis Matrix for shifting
       open(unit=11,file = hfile, status = 'old')
       
       close(11)
! Read in the vector
       open(unit=12, file = vecfile, status = 'unknown')
       do 10 i = 1, veclength
          read(10, *) vec(i)
10     end do
       close(12)
       do 20 i= row_start,row_stop
       
20     end do

       end
