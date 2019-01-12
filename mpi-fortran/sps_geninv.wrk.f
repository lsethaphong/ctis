!======================================================================
! FILENAME:
!    sps_geninv.f
!
! DESCRIPTION:
!    This code gives the generalized inverse of the hbasis projection
!    matrix stored as a series of pixel indexed files in signed sparse 
!    row format.
!    Adapted for very large sparse matrix in Fortran-MPI from
!    Fast Computation of Moore-Penrose Inverse Matrices
!    by Pierre Courrieu
!    Neural Information Processing -- Letters and Reviews
!    Vol. 8, No.2, August 2005
!
! HISTORY:
!   November 17, 2005: Original Code
!   April 5, 2006    : Completed generation of matrix A
!                      sps_HtrxH
!   April 9, 2006    : Completed generation of matrix L
!
!   April 12, 2006   : Inversion of L
!  
!   April 13, 2006   : Y=L*M*M*L'*G' multiplication routines
!                      where M = inv(L'*L)
!
!   Author:           Latsavongsakda Sethaphong
!   Institution:      Molecular Physiology & Biophysics
!                     Piston Lab
!                     Vanderbilt University
!                     Nashville, TN
!                     latsavongsakda.sethapong@vanderbilt.edu
!----------------------------------------------------------------------
       include 'sps_COM.f'
       program sps_geninv
       use com_methods
       IMPLICIT NONE
       include 'mpif.h'
! Basis matrix projection and inverse
       REAL, DIMENSION(h_n) :: diagA
       REAL, DIMENSION(h_n) :: diagL
       integer status(MPI_STATUS_SIZE), ierr
       REAL, DIMENSION(max_hbasis,2) :: hbasis
       LOGICAL sps_stop, sps_trans 
       INTEGER i,j,k,taskid,numtasks,sz_hbasis,sz_hpinv,numworkers
       INTEGER m,n,x_dim,y_dim,vxl_depth,r
       real e,etime,t(2), ltol
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
!      function Y = geninv(G)
!% Returns the Moore-Penrose inverse of the argument
!% Transpose if m < n
!      follow the rules of vxl_depth
! 
!  read in basis matrix hbasis.dat

       open(unit=10, file='hbasis.dat',status='old')
       i = 1
       do while(.TRUE.)
          read(10,end=45,fmt=*) hbasis(i,1), hbasis(i,2)
          i=i+1
       enddo
45     close(10)

!       open(unit=11, file='parameters.dat',status='old')
!          read(11,end=50,fmt=*) x_dim, y_dim, vxl_depth  
!50     close(11)
       sps_trans =.false.
!=====================================================================
! BEGIN MAIN TASKID SPECIFIC FUNCTIONS ###############################
!---------------------------------------------------------------------
!       gonow = .false.

       if (taskid.eq.MASTER) then
          print *, 'READ in HBASIS'
          e=etime(t)
       endif 

       numworkers = numtasks - 1
!=============================================
!  CONSTRUCT A = L'L
!---------------------------------------------

       open(unit=11, file='diagA.dat',status='old')
       i = 1
       do while(.TRUE.)
          read(11,end=65,fmt=*) diagA(i)
          i=i+1
       enddo
65     close(11)

!       CALL sps_HtrxH(hbasis,taskid,numworkers,diagA)
       ltol = MINVAL(diagA)*1e-9 ! new tolerance

!=============================================
!  CALL BLOCKING ROUTINE
!---------------------------------------------
       CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

!=============================================
!  PERFORM FULL CHOLESKY DECOMPOSITION
!---------------------------------------------

!       CALL sps_CHOLC(taskid,numworkers,ltol)

!=============================================
!  PERFORM inverse(L)
!---------------------------------------------
       open(unit=11, file='xdata/diagL.dat',status='old')
       i=1
       do while(.TRUE.)
          read(11,end=75,fmt=*) diagL(i)
          i=i+1
       enddo
75     close(11)

       CALL sps_Linv(diagL,taskid,numworkers,ltol)

!=============================================
!  CALL BLOCKING ROUTINE
!---------------------------------------------
       CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

!=============================================
!  PERFORM MULTPILCATION TO GET Y
!  WHERE M = Linv * Linv' = inv(L'L)
!---------------------------------------------
!       if (sps_trans.eq.true) then
!          Y=G'*L*M*M*L';
!       else
!          Y=L*M*M*L'*G';
!       endif
!=============================================

!       CALL sps_Y(taskid,numworkers,hbasis,fnameY,sps_trans)

!=============================================
! END MAIN PROGRAM
!---------------------------------------------
       if (taskid.eq.MASTER)then
          print *, 'elapsed:', e, 'user:', t(1), 'sys:', t(2)
          print *, 'COMPLETED PROGRAM'
       endif
200    end 

!=============================================
! END ########################################
!---------------------------------------------

!=====================================================================
! ROUTINE NAME: sps_HtrxH
!
! DESCRIPTION:
!   This routine returns the matrix A = H'*H to be used in the Cholesky
!   full factorization routine to generate the lower triangular L.
!   We assume that m > n for an H_mn matrix. If m < n, then A = H*H'.
!   Adapted from the backward projection routine of SPS_RECON since
!   that is an n*m routine.
!   Construct the matrix A, one column at a time.  Only the lower
!   triangle of A is needed for the full decomposition. Hence, only
!   this part will be stored by column in row sparse format
!
! HISTORY:
!   December 7, 2005  : Adapted for MPI in large format
!
!   Author:           Latsavongsakda Sethaphong
!   Institution:      Molecular Physiology & Biophysics
!                     Piston Lab
!                     Vanderbilt University
!                     Nashville, TN
!                     latsavongsakda.sethapong@vanderbilt.edu
!---------------------------------------------------------------------

       subroutine sps_HtrxH(hbasis,taskid,numworkers,diagA)
       use com_methods
       IMPLICIT NONE
       include 'mpif.h'
!       real, dimension(h_n) :: f_k
!       real, dimension(h_n) :: g_k
       real, dimension(max_hbasis,2), intent(IN) :: hbasis
       integer, intent(IN) :: taskid
       integer, intent(IN) :: numworkers
!       real, intent(out) :: tol
       real, dimension(h_n), intent(out) :: diagA
       real, dimension(h_n) :: tempD,tempM
       INTEGER STATUS(MPI_STATUS_SIZE)
       integer tmp,i,j,summ,mtype, ierr
       integer srt_pxl,stp_pxl,yline_pxl,extra

! allocate by pixel lines of length src_x
!determine extra pixels
!       extra = mod((numworkers+1),h_n/vxl_cnt)
!       avg_pxl = (h_n/vxl_cnt-extra)/(numworkers+1)
        extra = mod((numworkers+1),src_y) ! extra y lines not covered
        yline_pxl = (src_y-extra)/(numworkers + 1) 
       if (extra.gt.0.and.numworkers.gt.0) then
          extra = extra + 1
       endif

       diagA =0.

!taskid starts at zero
       if (taskid.eq.MASTER) then
!          PRINT *, 'called sps_HtrxH'
          srt_pxl = 1
          stp_pxl = ((taskid+1)*yline_pxl + extra)*src_x
       else
          srt_pxl = (taskid*yline_pxl + extra)*src_x + 1
          stp_pxl = ((taskid+1)*yline_pxl + extra)*src_x
       endif
!          PRINT *, 'task', taskid, 'srt_pxl',srt_pxl,'stp_pxl',stp_pxl

!       PRINT *, 'made it here'
!======================================================================
       call CORE_HtrxH(srt_pxl,stp_pxl,hbasis,tempD)
!======================================================================
!   send completion message back to main
!======================================================================
       mtype = PCI_AK_R
!       print *, 'finuto '
!       return
       if (taskid.eq.MASTER) then
!          tmp = 0
!          mtype = PCI_AK_R
          do j = 1,numworkers
!            call MPI_RECV(tmp,1,MPI_INTEGER, j, mtype,
!     .             MPI_COMM_WORLD, status, ierr)
            call MPI_RECV(tempM,h_n,MPI_REAL, j, mtype,
     .             MPI_COMM_WORLD, status, ierr)
             tempD=tempM + tempD
          enddo
       print *, 'constructed A'
           ! save the diagonal file
           ! doesn't have to be in the main data folders
       open(unit=14, file='diagA.dat',status='replace',
     ~      position='append')
          do j = 1, h_n
             write(14,*) tempD(j)
          enddo
       close(14)
           diagA=tempD
       else ! waiting protocol
!        tmp=1
!        call MPI_SEND(tmp,1, MPI_INTEGER, MASTER,
!     .       mtype, MPI_COMM_WORLD, ierr)
        call MPI_SEND(tempD,h_n,MPI_REAL, MASTER,
     .       mtype, MPI_COMM_WORLD, ierr)
       endif

       call MPI_BCAST(diagA,h_n, MPI_REAL,
     .    MASTER, MPI_COMM_WORLD, ierr)
!===============================================
!  END OF H'H ##################################
!-----------------------------------------------
       end subroutine sps_HtrxH


!======================================================================
!  ROUTINE: CORE_HtrxH
!  
!  DESCRIPTION:  Support routine of distributed generation of H'H
!
! HISTORY:
!   December 7, 2005: Original Code
!
!   Author:           Latsavongsakda Sethaphong
!   Institution:      Molecular Physiology & Biophysics
!                     Piston Lab
!                     Vanderbilt University
!                     Nashville, TN
!                     latsavongsakda.sethapong@vanderbilt.edu!  
!----------------------------------------------------------------------
       subroutine CORE_HtrxH(srt_pxl,stp_pxl,hbasis,diagA)
       use com_methods
! file is saved by col_id
! need to get pixel id and corresponding row id of f_p
       IMPLICIT NONE
!       real, dimension(h_n), intent(OUT) :: h_i
       real, dimension(h_m):: h_j
       real, dimension(h_n):: a_ij
       integer, intent(IN) :: srt_pxl
       integer, intent(IN) :: stp_pxl
       real, dimension(max_hbasis,2), intent(IN) :: hbasis
       real, dimension(h_n), intent(OUT) :: diagA
       real, dimension(max_hbasis,2) :: hbasis_i !shifted i
       real, dimension(max_hbasis,2) :: hbasis_j ! shifted j 
       integer, dimension(vxl_cnt+1) :: hrowi, hrowj, hrow
       real summ
       integer k,l,m,n,x_pos,y_pos,j,i_l,i_m,i_k,i
       character (LEN=24) :: fnameA, numAsString

        call strip(hrow,hbasis,vxl_cnt) ! gets shifted row/col
        do 200 k = srt_pxl, stp_pxl ! which pixels between src_x, src_y
           x_pos = MOD(k,src_x) ! src_x is defined in com_methods
           y_pos = (k-x_pos)/src_x + 1 ! get y position

           if (x_pos.eq.0) then ! final check on x position
              x_pos = src_x
              y_pos = y_pos-1
           endif

           write(numAsString,'(i0)') y_pos
           fnameA = 'xdata/a'//TRIM(numAsString)//'.dat' ! concat

           j = (k-1)*vxl_cnt ! gets starting voxel row

!           print *, j, ' column'
!           print *, 'x',x_pos,'y',y_pos
           call spsMAT_shift(hbasis,x_pos-1,y_pos-1,fpa_x,hbasis_j)
 
! get h_j loaded           

           do 220 m = 1, vxl_cnt ! decompse to voxels
              h_j = 0. ! clear h_j entries
              a_ij = 0. ! start new j column of a_ij
              summ =0.
              do 210 l = (hrow(m)+1), hrow(m+1) ! indices of hbasis j
                 h_j(INT(ABS(hbasis_j(l,2)))) = hbasis_j(l,1)
                 summ = hbasis_j(l,1)+summ
210           continue
!              print *, 'magnitude of hbasis_j',summ
!              return ! break
!            get h_jj which is the diagonal term  of a_ii
              do 235 i_m = m, vxl_cnt
                 summ = 0. ! clear for next sum i_m
                 do 240 i_l = (hrow(i_m)+1),hrow(i_m+1)
                    summ = summ +
     ~       hbasis_j(i_l,1)*h_j(INT(ABS(hbasis_j(i_l,2))))
240              continue !i_l
                 if (summ.gt.tol) then
                    a_ij(j+i_m) = summ
                 endif
235           continue ! finished first a_ij = h_i'*h_j

!====================================================================
! end of H_i=H_j
! now loop over all h_i's where i <= j, starting with the same pixel
!====================================================================
! do the rest of the pixels saving only the lower triangular
!====================================================================
               do 225 i_k = k+1, src_x*src_y ! go to the end of column
                  x_pos = MOD(i_k,src_x)
                  y_pos = (i_k-x_pos)/src_x + 1
                  if (x_pos.eq.0) then
                    x_pos = src_x
                    y_pos = y_pos-1
                  endif
                  n = (i_k-1)*vxl_cnt ! update new voxel start
               call spsMAT_shift(hbasis,x_pos-1,y_pos-1,fpa_x,hbasis_i)
                  do 245 i_m = 1, vxl_cnt
                     summ = 0. ! clear for next sum i_m
                     do 250 i_l = (hrow(i_m)+1),hrow(i_m+1)
                        summ = summ +
     ~        hbasis_i(i_l,1)*h_j(INT(ABS(hbasis_i(i_l,2))))
250                  continue !i_l 
                     if (summ.gt.tol) then
                        a_ij(n+i_m) = summ
                     endif
245               continue !i_m row m

225            continue !i_k voxel group k

           open(unit=14, file=fnameA,status='old',ERR=201,
     ~     position='append')
           goto 202
201        open(unit=14, file=fnameA,status='new',position='append')
202          do 255 i = j+m, h_n ! full length of the column j in pixel k
                 if (i.eq.(j+m)) then
                   diagA(i)= a_ij(i)
                 endif
                 if (a_ij(i).ge.tol) then
                    write(14,*)  a_ij(i), i
                 endif
255           continue
              write(14,*) 0,-(j+m) ! voxel column j in order, delimiter
           close(14)
220        continue !m column for voxel group m of pixel k
        
200     continue
       end subroutine CORE_HtrxH
!======================================================================
!   END H'H routines ##################################################
!----------------------------------------------------------------------




!======================================================================
! ROUTINE NAME: sps_CHOLC
!
! DESCRIPTION:
!   By the method of P. Courrieu, this routine performs the full rank
!   Cholesky decomposition of a matrix A = LxL'. For this application, we
!   give the shift invariant basis, hbasis.  Output of L is in sparse
!   format and separate diagonal file diagL.dat
!   The dimnesion of L is the same as  A.
!       L(k:n,r) = A(k:n,k)-L(k:n,1:(r-1)*transpose(L(k,1:(r-1))
!           if (L(k,r).gt.tol) then
!              L(k,r) =sqrt(L(k,r))
!              if (k.lt.n)then ! lower bound of column" r" is thus
!                 L((k+1):n,r)=L((k+1):n,r)/L(k,r)
!              endif
!           else
!             r=r-1
!           endif
!
!   By numerical recipes Press et al.
!   SUBROUTINE choldc(a,n,np,p)
!   INTEGER n,np
!   REAL a(np,np),p(n)
!   Given a positive-definite symmetric matrix a(1:n,1:n), 
!   with physical dimension np, this
!   routine constructs its Cholesky decomposition, A = L·LT . 
!   On input, only the upper triangle
!   of a need be given; it is not modified. 
!   The Cholesky factor L is returned in the lower triangle
!   of a, except for its diagonal elements which are returned in p(1:n).
!   INTEGER i,j,k
!   REAL sum
!      do 13 i=1,n
!         do 12 j=i,n
!         sum=a(i,j)
!            do 11 k=i-1,1,-1
!               sum=sum-a(i,k)*a(j,k)
!            enddo 11
!            if(i.eq.j)then
!               if(sum.le.0.)pause 'choldc failed'!a, with rounding errors, 
!                                                 !is not positive definite.
!                  p(i)=sqrt(sum)
!               else
!                  a(j,i)=sum/p(i)
!               endif
!            enddo 12
!         enddo 13
!      return
!   END
!
! HISTORY:
!   December 7, 2005: Original Code
!
!   Author:           Latsavongsakda Sethaphong
!   Institution:      Molecular Physiology & Biophysics
!                     Piston Lab
!                     Vanderbilt University
!                     Nashville, TN
!                     latsavongsakda.sethapong@vanderbilt.edu
!----------------------------------------------------------------------

       subroutine sps_CHOLC(taskid,numworkers,ltol)
       use com_methods
       IMPLICIT NONE
       include 'mpif.h'
       INTEGER, INTENT(IN) :: taskid
       INTEGER, INTENT(IN) :: numworkers ! numtasks-1
       REAL, INTENT(IN) :: ltol
!       INTEGER, INTENT(IN) :: k
!      FORTRAN ARRAYS ARE NUMBERED 1 up....
!       REAL, DIMENSION(h_n) :: A_k ! voxel element value array
!       INTEGER, DIMENSION(h_n) :: A_idx ! voxel per x line
       REAL, DIMENSION(h_n) :: A_val ! assume max non zero is src_y
       REAL, DIMENSION(h_n) :: L_r  ! must be continually updated
       REAL, DIMENSION(h_n) :: L_s  ! subtractand vector
       REAL, DIMENSION(h_n) :: diagL ! diagonal of L
       REAL val
       INTEGER, DIMENSION(2) :: cols
       INTEGER STATUS(MPI_STATUS_SIZE)
       INTEGER  i_r,k,sk,r,ierr,mtype,i_j
       INTEGER i_k,i,it,ip,fk,j,v,i_s,idx
       CHARACTER (LEN=24) :: fnameA, numAsString

!==============================================
!   Load A by pixel voxel group and operate a 
!   a column of A at a time
!
!   The MASTER routine handles the passing
!   protocols and main loop
!----------------------------------------------
       if (taskid.eq.MASTER) then
          print *, 'Cholesky Factoring L'
          print *, 'tolerance is',ltol
          k   = 0
          r   = 0  ! initialize r selected column of L
          L_s = 0. ! initialize subtractand L
          diagL = 0. ! initialize diagonal array of L
       do 410 i_k = 1, src_y ! A_k files are grouped by y lines
!==============================================
!    update row file of Ak           
!----------------------------------------------
          write(numAsString,'(i0)') i_k
          fnameA = 'xdata/a'//TRIM(numAsString)//'.dat' ! concat
!          print *, 'opening',fnameA
          open(unit=14, file=fnameA,status='old')
!===============================================
!  STOP LOAD
!-----------------------------------------------
             print *, 'reading ',fnameA
             j = 1 ! initialize j
          do 415 v = 1,src_x*vxl_cnt ! A_idx thru a row of x voxels
!==============================================
!    update row file of Ak
!----------------------------------------------
             i     = 1
             A_val = 0.
             val = 0.
             idx = 0
             do i_j =1,h_n
                read(14,end=445,fmt=*,ERR=446) val,idx
                if (idx.lt.0) then
                   goto 409
                elseif (i.gt.h_n) then
                   print *, 'reading file error',fnameA
                   close(14)
                   return
                endif
                A_val(idx)=val
                i=i+1
             enddo
445       close(14)
          goto 409
446       print *, 'cannot read',fnameA   
!===============================================
!  STOP LOAD
!-----------------------------------------------
409           k=k+1   ! the column id
!             A_val = 0. ! clear A_val
!             do while(A_k(j).gt.0)
!                A_val(A_idx(j)) = A_k(j)
!                j=j+1 ! advance j
!             enddo
!             j=j+1 ! advance j
             r=r+1 ! advance r selected column of L
!===============================================
!  CONSTRUCT L(k:n,r) 1 column of r at a time
!----------------------------------------------- 
             L_r = 0. ! clear L_r
             L_r(k:h_n)=A_val(k:h_n)-L_s(k:h_n)
             if (L_r(k).gt.ltol) then
                L_r(k)=SQRT(L_r(k))
                diagL(k)=L_r(k) ! save the diagonal element 
                if(k.lt.h_n) then
                   L_r((k+1):h_n)=L_r((k+1):h_n)/L_r(k)
                endif
             else
                r=r-1 ! stay at this column for the next round
             endif
!            L(1:n,1) is A(1:n,1) first
!===============================================
!  WRITE TO FILES OF Lxxx.dat grouped by i_k
!-----------------------------------------------
           PRINT *, 'filing',i_k
           write(numAsString,'(i0)')i_k 
           fnameA='xdata/L'//TRIM(numAsString)//'.dat' ! Lfile
           open(unit=15,file=fnameA,status='old',ERR=401,
     ~     position='append')
           goto 403
401        open(unit=15,file=fnameA,status='new',ERR=402,
     ~     position='append')
           goto 403
!===============================================
402        print *, 'cant open',fnameA
           return
!-----------------------------------------------
403             do 420 i_r=k,h_n 
                   if (L_r(i_r).gt.ltol) then
                      write(15,*) L_r(i_r), i_r
                   endif
420             continue
               write(15,*) 0,-k
            close(15)
!===============================================
!   ############################################
!-----------------------------------------------
!===============================================
!  CONSTRUCT L_s(k+1:n,1:r)=L(k+1:n,1:r)*L(k+1,1:r)'
!  BASED ON THE NEXT EXPECTED r column of L
!-----------------------------------------------
            if (k.lt.h_n) then
               L_s=0. ! clear subtractand
               L_r=0. ! reuse old vector
!===============================================
!  CALL BARRIER BEFORE PROCEEDING
!-----------------------------------------------
!               CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
               cols=(/k+1,r/)
               mtype = PCI_LS_S
               do i_s = 1,numworkers
                   if (r.gt.i_s) then
               CALL MPI_SEND(cols,2,MPI_INTEGER,i_s,
     .             mtype, MPI_COMM_WORLD, ierr)
                   endif
               enddo
!               CALL MPI_BCAST(cols,2,MPI_INTEGER,MASTER,
!     .              MPI_COMM_WORLD, ierr)
               sk = cols(1)
               fk = h_n-sk+1
!===============================================
!  RETRIEVE L_s sums from slaves
!-----------------------------------------------
               CALL LS_PARTIAL_SUM(sk,r,L_s,taskid,numworkers)
               mtype = PCI_LS_R
               do i_s = 1,numworkers
                  if (r.gt.i_s) then
!                  print *, 'calling Ls from',i_s
               CALL MPI_RECV(L_r(sk:h_n),fk,MPI_REAL,i_s,
     .                 mtype,MPI_COMM_WORLD, status, ierr)
               L_s(sk:h_n)=L_s(sk:h_n)+L_r(sk:h_n)
                  endif
               enddo
            endif

415       continue  
!===============================================

!-----------------------------------------------
410    continue ! main i_k loop/MASTER LOOP
!===============================================
!  Save diagonal of L in diagL file
!-----------------------------------------------
          open(unit=17,file='xdata/diagL.dat',status='new',
     ~        position='append')
          do i=1,h_n
             write(17,*) diagL(i)
          enddo
          close(17)
!===============================================
!  WORKER ROUTINES #############################
!-----------------------------------------------
       else
!===============================================
!  ROUTINE to generate L_s
!-----------------------------------------------
!===============================================
!  WAIT FOR THE MASTER REQUEST
!-----------------------------------------------
              k=1
           do while(k.lt.h_n)
              L_s = 0.
              L_r = 0.
!===============================================
!   WAIT FOR THE MASTER REQUEST
!-----------------------------------------------
!              CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
!              CALL MPI_BCAST(cols,2,MPI_INTEGER,MASTER,
!     .             MPI_COMM_WORLD, ierr)
              mtype = PCI_LS_S
              CALL MPI_RECV(cols,2,MPI_INTEGER,MASTER,
     .             mtype,MPI_COMM_WORLD, status, ierr)
              k = cols(1) ! row id which is "k+1"
              r = cols(2) ! last column
              if (k.lt.1) then
                 print *, 'invalid k value exiting'
                 return
!              elseif(k.eq.taskid) then
!                 print *, taskid,'doing well'
              endif
              fk = h_n-k+1 ! number of rows to return
!              CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
              CALL LS_PARTIAL_SUM(k,r,L_r,taskid,numworkers)
              mtype = PCI_LS_R ! send back on the relevant rows
              CALL MPI_SEND(L_s(k:h_n),fk,MPI_REAL,MASTER,
     .             mtype, MPI_COMM_WORLD, ierr)
!              CALL MPI_REDUCE(L_r(k:h_n),L_s(k:h_n),fk,MPI_REAL,
!     .               MPI_SUM,MASTER,MPI_COMM_WORLD,ierr)
!              print *, taskid, 'sent Ls', k, 'for r',r             
           enddo ! Wait Loop
        endif
!        CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
!        print *, taskid, 'has exited CHOLC'
        return
!=====================================================================
!  END OF sps_CHOLC ##################################################
!---------------------------------------------------------------------
       end subroutine sps_CHOLC

!=====================================================================
! ROUTINE NAME: LS_PARTIAL_SUM
!
! DESCRIPTION:
!   This routine returns the partial sum of the L subtractand
!   in the generation of the lower triangular matrix in support
!   of the full rank Cholesky factorization adapated for multi
!   processor distribution and file handling.
!
! HISTORY:
!   April 9, 2006  : Adapted for MPI in large format
!
!   Author:           Latsavongsakda Sethaphong
!   Institution:      Molecular Physiology & Biophysics
!                     Piston Lab
!                     Vanderbilt University
!                     Nashville, TN
!                     latsavongsakda.sethapong@vanderbilt.edu
!---------------------------------------------------------------------
       subroutine LS_PARTIAL_SUM(k,r,L_s,taskid,numworkers)
       use com_methods
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: k ! row element to scale columns
       INTEGER, INTENT(IN) :: r ! 1 to r columns needing partial sums
       REAL, DIMENSION(h_n), INTENT(OUT) :: L_s ! summing vector
       INTEGER, INTENT(IN) :: taskid
       INTEGER, INTENT(IN) :: numworkers ! numtasks -1
       REAL, DIMENSION(h_n) :: L_k ! y voxel column storage array
       REAL, DIMENSION(h_n) :: L_c ! column to be scaled by element k
       INTEGER, DIMENSION(h_n) :: L_idx ! index storage array of L_k values
       INTEGER, DIMENSION(src_x*vxl_cnt) :: L_vxl !voxel index pointer array
       INTEGER i,j,idx,idc, y_file, srt_yfile, stp_yfile
       INTEGER extra, avg_col, r_pos, y_pos, i_vxl
       INTEGER srt_col, stp_col, i_c
       REAL val
       CHARACTER (LEN=24) :: fnameA, numAsString
!===============================================
!   L_r(fk:h_n)=L_r(k:h_n)+L_k(k)*L_k(k:h_n)
!-----------------------------------------------
!===============================================
!   SELECT RANGE OF Y FILES BASED ON TASKID 
!   AND EXPECTED NUMBER OF EXISTENT Y FILES
!   THEN PRE SELECT COLUMNS
!----------------------------------------------- 
!===============================================
!   split files based on taskid, numtasks, and r
!-----------------------------------------------
!     determine extra columns       
          extra = mod((numworkers+1),r)
          avg_col = (r-extra)/(numworkers+1)
          if (taskid.ge.r) then
             print *, taskid, 'no deal on',r
             L_s = 0.
             return
          endif
!      extra = r - avg_col*numtasks
!          print *, taskid, ' called partial sum for', r
!taskid starts at zero
          if (taskid.eq.MASTER) then
!             print *, 'Allocating columns',r,'avg ',avg_col
             srt_col = 1
             stp_col = (taskid+1)*avg_col + extra
!             print *, 'start',srt_col,'stop',stp_col
          else
             srt_col = taskid*avg_col + extra + 1
             stp_col = (taskid+1)*avg_col + extra
          endif
! determine first and last y files

           r_pos = MOD(srt_col,src_x*vxl_cnt) ! src_x is defined in com_methods
           y_pos = (srt_col-r_pos)/(src_x*vxl_cnt) + 1 ! get y position

           if (r_pos.eq.0) then ! final check on x position
!              r_pos = src_x*vxl_cnt ! last col of file y group
              y_pos = y_pos-1
           endif
!
          srt_yfile = y_pos

          r_pos = MOD(stp_col, src_x*vxl_cnt)
          y_pos = (stp_col-r_pos)/(src_x*vxl_cnt) + 1

          if (r_pos.eq.0) then
             y_pos = y_pos-1
!             print *, y_pos, '+++r_pos is zero', stp_col
          endif

          stp_yfile = y_pos
!===============================================
!  END of y_file determination, loop through
!  file data and perform partial summation
!-----------------------------------------------
          L_s = 0.
          do y_file = srt_yfile, stp_yfile
          i_vxl = (y_file-1)*src_x*vxl_cnt ! gives the major voxel index
          write(numAsString,'(i0)')y_file
          fnameA = 'xdata/L'//TRIM(numAsString)//'.dat' ! concat
!          print *, taskid,'opening ',fnameA
          open(unit=19,file=fnameA,status='old')
!===============================================
!  Open file read protocols
!-----------------------------------------------
455       L_c = 0.  ! real
          idx = 0
          val = 0.
          i = 1 ! just an index counter
          do while(.TRUE.)
             read(19,end=445,fmt=*) val, idx
             if (idx.lt.0) then
                i_c = ABS(idx)
                goto 450
             elseif (i.gt.h_n) then
                print *, 'file error'
                goto 445
                return
             endif
             L_c(idx)=val
             i=i+1
          enddo ! close the file
!445       print *, taskid,'closing ',fnameA
445       close(19)
!          print *, taskid,'closed ',fnameA
          goto 465
!===============================================
!  Ensure that the column id is within the selection
!-----------------------------------------------
450       if ((i_c.ge.srt_col).and.(i_c.le.stp_col)) then
                L_s = L_s + L_c(k)*L_c
              goto 455
          endif 
465    enddo ! get next y_file data
!       print *,taskid,'exited partial sum'
       return ! end of routine
       end subroutine LS_PARTIAL_SUM

!===============================================
!  END LS_PARTIAL_SUM by k element of column 
!-----------------------------------------------
!======================================================================
!  END OF CHOLC #######################################################
!----------------------------------------------------------------------





!======================================================================
! ROUTINE NAME: sps_Linv
!
! DESCRIPTION:
!   Take the outputs from cholc and allow for direct inversion of L
!   From Numerical Recipes in C++
!   Linv is:
!       for (i=0;i<n;i++) {
!          a[i][i]=1.0/p[i]; //p[i] diagonal elements of a
!          for(j=i+1;j<n;j++) {
!             sum=0.0;
!             for (k=i;k<j;k++) sum -= a[j][k]*a[k][i];
!             a[j][i]=sum/p[j];
!          }
!       }
!   From Numerical Recipes in Fortran
!   Linv is 
!      do 13 i=1,n
!         a(i,i)=1./p(i)
!         do 12 j=i+1,n
!            sum=0.
!            do 11 k=i,j-1
!               sum=sum-a(j,k)*a(k,i)
!            enddo 11
!            a(j,i)=sum/p(j)
!         enddo 12
!      enddo 13
!   
!
! HISTORY:
!   January 11, 2006: Original Code
!
!   April 12, 2006  : Complete y_file inversion a column at a 
!                     time
!
!   Author:           Latsavongsakda Sethaphong
!   Institution:      Molecular Physiology & Biophysics
!                     Piston Lab
!                     Vanderbilt University
!                     Nashville, TN
!                     latsavongsakda.sethapong@vanderbilt.edu
!----------------------------------------------------------------------
       subroutine sps_Linv(p,taskid,numworkers,ltol)
       use com_methods
       IMPLICIT NONE
       include 'mpif.h'
       REAL, DIMENSION(h_n), INTENT(IN) :: p
       INTEGER, INTENT(IN) :: taskid
       INTEGER, INTENT(IN) :: numworkers ! numtasks-1
       REAL, INTENT(IN) :: ltol
!      FORTRAN ARRAYS ARE NUMBERED 1 up....
       REAL, DIMENSION(h_n) :: Linv ! assume max non zero is src_y
       REAL, DIMENSION(h_n) :: L_r  ! row vector of j
       REAL, DIMENSION(h_n) :: L_c  !  
       REAL val, summ
       INTEGER, DIMENSION(2) :: cols
       INTEGER STATUS(MPI_STATUS_SIZE)
       INTEGER  i_r,k,sk,r,ierr,mtype,i_j
       INTEGER i_k,i,it,ip,fk,j,v,i_s,idx
       CHARACTER (LEN=24) :: fnameA, numAsString

!==============================================
!   Load A by pixel voxel group and operate a
!   a column of A at a time
!
!   The MASTER routine handles the passing
!   protocols and main loop
!----------------------------------------------
       if (taskid.eq.MASTER) then
          print *, 'Inverting Cholescky factor L'
          print *, 'tolerance is',ltol
          k   = 0
       do 510 i_k = 1, src_y ! L files are grouped by y lines
!==============================================
!    read in successive columns of L
!----------------------------------------------
          write(numAsString,'(i0)')i_k
          fnameA = 'xdata/L'//TRIM(numAsString)//'.dat' ! concat
          open(unit=14, file=fnameA,status='old')
!===============================================
!  STOP LOAD
!-----------------------------------------------
!             print *, 'reading ',fnameA
!             j = 1 ! initialize j
!==============================================
!  Go through the columns in a file of y
!----------------------------------------------
          do 515 v = 1,src_x*vxl_cnt ! goes through columns of y
                i = 1  ! counting index
             Linv = 0. ! new column i 
              val = 0. ! temp value holder
              idx = 0  ! temp index holder
             do i_j =1,h_n
                read(14,end=545,fmt=*,ERR=546) val,idx
                if (idx.lt.0) then
                   k = ABS(idx) ! L column just loaded
                   goto 509 ! finished reading in column go to routine
                elseif (i.gt.h_n) then
                   print *, 'reading file error',fnameA
                   close(14)
                   return
                endif
                Linv(idx)=val
                i=i+1
             enddo
545       close(14) ! end of file
          goto 509 ! open next file
546          print *, 'cannot read',fnameA
!===============================================
!  STOP LOAD AND MAKE THINGS HAPPEN
!-----------------------------------------------
!===============================================
!  CONSTRUCT Linv k column at a time
!-----------------------------------------------
509       Linv(k)=1./p(k) ! diagonal term
          if (k.gt.h_n) then ! k must be > h_n
!===============================================
!  THE WORKERS WILL WAIT ON COLS
!  THE KILL SIGNAL IS k = - taskid
!-----------------------------------------------
          do j = k+1,h_n ! rows of j are  k+1 to the end
             cols=(/k,j/) ! send column and row
             do i_s = 1,numworkers
                if ((j-k).gt.i_s) then
                   mtype = PCI_LI_S
                   CALL MPI_SEND(cols,2,MPI_INTEGER,i_s,
     .                 mtype,MPI_COMM_WORLD,ierr)
                   mtype = PCI_LK_S
                   CALL MPI_SEND(Linv,h_n,MPI_REAL,i_s,
     .                 mtype,MPI_COMM_WORLD,ierr)
                endif
             enddo
!===============================================
!  GET PARTIAL SUM FROM MASTER FOR SPAN j - k
!-----------------------------------------------
             CALL CORE_Linv(k,j,Linv,summ,taskid,numworkers)
!===============================================
!  GET PARTIAL SUMS FROM SLAVES
!-----------------------------------------------
             mtype = PCI_LI_R
             do i_s = 1,numworkers
                ! call slaves only if sufficient span
                if ((j-k).gt.i_s) then
                   CALL MPI_RECV(summ,1,MPI_REAL,i_s,
     .                mtype,MPI_COMM_WORLD, status, ierr)
                  Linv(j)=Linv(j)+summ ! summ should be negative
                endif
             enddo ! end of row element j
             Linv(j)=Linv(j)/p(j)
          enddo ! of column k

!===============================================
!  k=h_n, end of inversion, stop slaves
!-----------------------------------------------
          else
             do i_s =1,numworkers
                   cols=(/-1,-1/) ! kill signal
                   mtype = PCI_LI_S
                   CALL MPI_SEND(cols,2,MPI_INTEGER,i_s,
     .                 mtype,MPI_COMM_WORLD,ierr)
             enddo
          endif ! k must be > h_n
!===============================================
!  END OF CONSTRUCTING COLUMN k of Linv
!-----------------------------------------------
 
!===============================================
!  APPEND COLUMN OF k into Linv files grouped by i_k
!-----------------------------------------------
          PRINT *, 'filing',i_k
          write(numAsString,'(i0)')i_k
          fnameA='xdata/Linv'//TRIM(numAsString)//'.dat' ! Lfile
          open(unit=15,file=fnameA,status='old',ERR=501,
     ~    position='append')
          goto 503
501       open(unit=15,file=fnameA,status='new',ERR=502,
     ~    position='append')
          goto 503
!===============================================
502       print *, 'cant open',fnameA
          return
!-----------------------------------------------
503          do 520 i_r=k,h_n
                if (Linv(i_r).gt.ltol) then
                   write(15,*) Linv(i_r), i_r
                endif
520          continue
             write(15,*) 0,-k
             close(15)
!===============================================
!-----------------------------------------------
515       continue ! voxels per y_file
!-----------------------------------------------
510    continue ! main i_k loop/MASTER LOOP
!===============================================
! END OF MASTER ROUTINE FOR Linv
!-----------------------------------------------

!===============================================
! **********************************************
!  WORKER ROUTINES #############################
! **********************************************
!-----------------------------------------------
       else
!===============================================
!  ROUTINE to generate L_s
!-----------------------------------------------
!===============================================
!  WAIT FOR THE MASTER REQUEST
!-----------------------------------------------
              k=1
           do while(k.lt.h_n)
!===============================================
!   WAIT FOR THE MASTER REQUEST
!-----------------------------------------------
              Linv = 0. ! clear Linv
              mtype = PCI_LI_S
              CALL MPI_RECV(cols,2,MPI_INTEGER,MASTER,
     .             mtype,MPI_COMM_WORLD, status, ierr)
              k = cols(1) ! current column of Linv
              j = cols(2) ! row element
              if (k.lt.1) then
                 print *, taskid,'exiting Linv'
                 return
              endif
              mtype = PCI_LK_S
              CALL MPI_RECV(Linv,h_n,MPI_REAL,MASTER,
     .             mtype,MPI_COMM_WORLD, status, ierr)
!===============================================
!  Get the partial sum from the general routine
!-----------------------------------------------
              CALL CORE_Linv(k,j,Linv,summ,taskid,numworkers)
!===============================================
!  SEND BACK PARTIAL SUM
!-----------------------------------------------
              mtype = PCI_LI_R ! send back on the relevant rows
              CALL MPI_SEND(summ,1,MPI_REAL,MASTER,
     .             mtype, MPI_COMM_WORLD, ierr)
           enddo ! Wait Loop
        endif
        return ! exit routine
!=====================================================================
!  END OF sps_Linv  ##################################################
!---------------------------------------------------------------------
       end subroutine sps_Linv

!======================================================================
! ROUTINE NAME: CORE_Linv
!
! DESCRIPTION:
!   From Numerical Recipes in Fortran
!   Linv is
!      do 13 i=1,n
!         a(i,i)=1./p(i)
!         do 12 j=i+1,n
!            sum=0.
!            do 11 k=i,j-1
!               sum=sum-a(j,k)*a(k,i)
!            enddo 11
!            a(j,i)=sum/p(j)
!         enddo 12
!      enddo 13
!   
!   This routine is performs the growing row sums for each elent Linv
!   to the task routine and saves the resultant column.  Files are
!   organized by y_line pixels to prevent common write.
!   The routine is adapted from LS_PARTIAL_SUM of sps_CHOLC
!
! HISTORY:
!   January 11, 2006: Original Code
!
!   Author:           Latsavongsakda Sethaphong
!   Institution:      Molecular Physiology & Biophysics
!                     Piston Lab
!                     Vanderbilt University
!                     Nashville, TN
!                     latsavongsakda.sethapong@vanderbilt.edu
!----------------------------------------------------------------------
       subroutine CORE_Linv(k,r,Linv,summ,taskid,numworkers)
       use com_methods
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: k ! current column
       INTEGER, INTENT(IN) :: r ! row for calculation
       REAL, DIMENSION(h_n), INTENT(IN) :: Linv ! growing inverse column
       REAL, INTENT(OUT) :: summ ! partial dot product
       INTEGER, INTENT(IN) :: taskid
       INTEGER, INTENT(IN) :: numworkers ! numtasks - 1
       INTEGER, DIMENSION(h_n) :: L_c ! row element of columns
       INTEGER i,j,idx,idc, y_file, srt_yfile, stp_yfile
       INTEGER extra, avg_col, r_pos, y_pos, i_vxl
       INTEGER srt_col, stp_col,i_c, span
       REAL val
       CHARACTER (LEN=24) :: fnameA, numAsString
!===============================================
!   COLUMN BY COLUMN GENERATION OF L inverse
!   INDIVIDUALIZED BY COLUMN AND TASK ID
!-----------------------------------------------
!===============================================
!   SELECT RANGE OF Y FILES BASED ON TASKID
!   AND EXPECTED NUMBER OF EXISTENT Y FILES
!   THEN PRE SELECT COLUMNS
!-----------------------------------------------
!===============================================
!   split files based on taskid, numtasks, k, and r
!-----------------------------------------------
          span = r-k 
!     determine extra columns
          extra = mod((numworkers+1),span) ! average number
          avg_col = (span-extra)/(numworkers+1)
          summ = 0.
          if (taskid.ge.span) then
             print *, taskid, 'no deal on',r
             return
          endif
!      extra = r - avg_col*numtasks
!          print *, taskid, ' called partial sum for', r
!taskid starts at zero
          if (taskid.eq.MASTER) then
             print *, 'Linv span ',span,'avg ',avg_col
             srt_col = k !starting column
             stp_col = k+(taskid+1)*avg_col + extra
!             print *, 'start',srt_col,'stop',stp_col
          else
             srt_col = k+taskid*avg_col + extra + 1
             stp_col = k+(taskid+1)*avg_col + extra
          endif
! determine first and last y files

           r_pos = MOD(srt_col,src_x*vxl_cnt) ! src_x is defined in com_methods
          if (taskid.eq.MASTER) then
!             print *, 'Allocating columns',r,'avg ',avg_col
             srt_col = k
             stp_col = k+(taskid+1)*avg_col + extra
!             print *, 'start',srt_col,'stop',stp_col
          else
             srt_col = k+taskid*avg_col + extra + 1
             stp_col = k+(taskid+1)*avg_col + extra
          endif
! determine first and last y files

           r_pos = MOD(srt_col,src_x*vxl_cnt) ! src_x is defined in com_methods
           y_pos = (srt_col-r_pos)/(src_x*vxl_cnt) + 1 ! get y position

           if (r_pos.eq.0) then ! final check on x position
!              r_pos = src_x*vxl_cnt ! last col of file y group
              y_pos = y_pos-1
           endif
!
          srt_yfile = y_pos

          r_pos = MOD(stp_col, src_x*vxl_cnt)
          y_pos = (stp_col-r_pos)/(src_x*vxl_cnt) + 1

          if (r_pos.eq.0) then
             y_pos = y_pos-1
!             print *, y_pos, '+++r_pos is zero', stp_col
          endif

          stp_yfile = y_pos
!===============================================
!  END of y_file determination, loop through
!  file data and perform partial summation
!-----------------------------------------------
          do y_file = srt_yfile, stp_yfile
          i_vxl = (y_file-1)*src_x*vxl_cnt ! gives the major voxel index
          write(numAsString,'(i0)')y_file
          fnameA = 'xdata/L'//TRIM(numAsString)//'.dat' ! concat
!          print *, taskid,'opening ',fnameA
          open(unit=19,file=fnameA,status='old')
!===============================================
655       L_c = 0. ! real column of L
          idx = 0  ! temp index holder
          val = 0. ! temp value holder
          i = 1 ! just a index counter
          do while(.TRUE.)
             read(19,end=645,fmt=*) val, idx
             if (idx.lt.0) then
                i_c = ABS(idx)
                goto 650
             elseif (i.gt.h_n) then
                print *, 'file error'
                goto 645
                return
             endif
             L_c(idx)=val
             i=i+1
          enddo ! close the file
645       close(19)
          goto 665
!===============================================
!  Ensure that the column id is within the selection
!-----------------------------------------------
650       if ((i_c.ge.srt_col).and.(i_c.le.stp_col)) then
                summ = summ - L_c(r)*Linv(i_c)
              goto 655
          endif
665    enddo ! get next y_file data
!       print *,taskid,'exited partial sum'
       return ! end of routine
       end subroutine CORE_Linv

!===============================================
!  END CORE_Linv by k element of column
!======================================================================
!  END OF Linv routiones###############################################
!----------------------------------------------------------------------

!======================================================================
! ROUTINE NAME: sps_Y
!
! DESCRIPTION:
!       Final iterative multiplication by the backward-forward
!       algorithm to give the generalized inverse of G
!       
!=============================================
!  PERFORM MULTPILCATION TO GET Y
!  WHERE M = Linv * Linv' = inv(L'L)
!---------------------------------------------
!       if (sps_trans.eq.true) then
!          Y=G'*L*M*M*L';
!       else
!          Y=L*M*M*L'*G'; <-- will calling this routine
!       endif
!
!       Y is saved as row successive format column sparse signed index
!
! HISTORY:
!   April 12, 2006:   Original Code
!
!   Author:           Latsavongsakda Sethaphong
!   Institution:      Molecular Physiology & Biophysics
!                     Piston Lab
!                     Vanderbilt University
!                     Nashville, TN
!                     latsavongsakda.sethapong@vanderbilt.edu
!----------------------------------------------------------------------
!       subroutine sps_Y(taskid,numtasks,tol)
!       INTEGER, INTENT(IN) :: taskid
!       INTEGER, INTENT(IN) :: numtasks
!       REAL, INTENT(IN) :: tol

!       end subroutine sps_Y

!       subroutine CORE_Y(taskid,numtasks,tol)
!       INTEGER, INTENT(IN) :: taskid
!       INTEGER, INTENT(IN) :: numtasks
!       REAL, INTENT(IN) :: tol

!       end subroutine CORE_Y
!======================================================================
!  END OF FILE ########################################################
!----------------------------------------------------------------------
