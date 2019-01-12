!================================================
! Image Reconstruction Main program
!================================================
! Upon Invocation, this ancillary program resident
! on the VAMPIRE system will read in the current
! h+ input file fed from the remote windows program.
!
! 
!
!
!
	program RECON()
!
! if check flag is high
! check for a base inversion h+ file that is up to date
! otherwise use the latest one

! call function to read in g -> g_o

! call function to perform initial guess of f from H+

! call expectation maximization function of f with g_o
! take the geometric difference

! call function to return f_best and send file vi sftp to 

! the requesting windows microscope operating program

	end

