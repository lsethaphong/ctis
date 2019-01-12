!======================================================================
! FILE NAME : sps_greville.f
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


!================================================
! from Jiri Rohn
! Institute of Computer Science
! Academy of Sciences of the Czech Republic
! Prehled nekterych dulezitych vet z teorie matic
! 5 Sep 2003
! Technical Report No. 893
! -------------------------
! modified for sparse matrix 
! The inverse will be saved by rowId to expedite
! multiplication
! -------------------------
! Utilizing Greville with Sparse matrices and vectors
! implementation of Greville
! 29 OCT 04
! Latsavongsakda Sethaphong
! Piston Lab
! VUMC
! -------------------------
! the system matrix is row sparse, hence, it is indexed by column
! hence, the inverse is column sparse, and, it is indexed by row.
! -------------------------
!===============================================

!===============================================
!  This is the basic Greville Routine
!
!
!===============================================
        subroutine GREVILLE(A, mi, ni, Hp)
        integer :: i,ii,iii,j,k,mi,ni
        dimension d(mi,1), c(mi,1), tmp_d(1,1),val_d(1,1), v_one(1,1)
        dimension X(ni,mi),dt(1,ni), bt(1,mi),tmp_X(ni,mi)
        dimension A(mi,ni)
        dimension Hp(ni,mi)
        real max_d,min_d,max_c,min_c,tol,val_d,v_one
        tol = 1.0e-6
        v_one = 1.0
        d(1:mi,1)=A(1:mi,1)! first column
        max_d = MAXVAL(d(1:mi,1))
        min_d = ABS(MINVAL(d(1:mi,1)))
        X=0.0
        if  (max_d.le.tol.or.min_d.le.tol) then
           bt=tol ! zeros
        else
           val_d = 1/(MATMUL(TRANSPOSE(d),d))
           bt=TRANSPOSE(MATMUL(d,val_d))
        endif
        X(1,1:mi)=bt(1,1:mi)
        do 101 j=2,ni
           d(1:(j-1),1)=MATMUL(X(1:(j-1),1:mi), A(1:mi,j))
           c(1:mi,1)=A(1:mi,j)-MATMUL(A(1:mi,1:(j-1)),d(1:(j-1),1))
           max_c = MAXVAL(c)
           min_c = MINVAL(ABS(c))
           if (max_c.lt.tol.or.min_c.lt.tol) then
                 val_d = 1.0
               do 102 k = 1,j-1
                 val_d  = val_d + d(k,1)**2
102            enddo
               val_d = 1/val_d
               dt=TRANSPOSE(d)! take only to the j-1th position
               bt(1,1:mi)=MATMUL(dt(1,1:(j-1)),X(1:(j-1),1:mi))
               bt=MATMUL(val_d,bt)
           else
              val_d = MATMUL(TRANSPOSE(c),c)
              val_d = 1/val_d
              bt=MATMUL(val_d,TRANSPOSE(c))
           endif
           tmp_X(1:j-1,:)=MATMUL(d(1:j-1,:),bt(:,:))
           X(1:j-1,1:mi)=X(1:j-1,1:mi)-tmp_X(1:j-1,1:mi)
           X(j:j,1:mi)=bt(:,:)
101      enddo
        Hp = X
        return
        end       
!===============================================
!
! Shift Expansion of Greville Inversion Theorem
! for sparse matrix calculation to get base 
! matrices, H, and Hinv in sparse row, col index
! format, respectively
!===============================================

!======================================================================
! ROUTINE NAME: sps_greville
!
! DESCRIPTION:
!
!
!
! HISTORY:
!   September 21, 2005: Original Code
!
!   Author:           Latsavongsakda Sethaphong
!   Institution:      Molecular Physiology & Biophysics
!                     Piston Lab
!                     Vanderbilt University
!                     Nashville, TN
!                     latsavongsakda.sethapong@vanderbilt.edu
!----------------------------------------------------------------------

       subroutine sps_greville(n_hcol,tol,ftol) 
       integer :: i,ii,iii,j,k,mi,ni
       dimension d(mi,1), c(mi,1), tmp_d(1,1),val_d(1,1), v_one(1,1)
       dimension X(ni,mi),dt(1,ni), bt(1,mi),tmp_X(ni,mi)
       dimension A(mi,ni)
       dimension Hp(ni,mi)
       real max_d,min_d,max_c,min_c,tol,val_d,v_one
       
       char*10 dirname
       char*60 fnameA, fnameB, fnameC, fnameD, fnameX, fnameK

! spsMat(colId,rowId,fpa_x,fpa_y,src_x,src_y,vxl_depth,nz_pts); % column format
! prod = sps_matmul(fnameA,fnameB,fnameC,func,b_ki,b_kn) % multiplication
! first calculate hbasis inverse
! then compute the corresponding shifted matrix in sparse format
! then use the partition theorem to grow the inverse matrix by 

       dirname = 'xdata/'
       fnameA = dirname + 'av'+ '.mat'
       fnameB = dirname + 'bv'+ '.mat'
       fnameC = dirname + 'cv'+ '.mat'
       fnameD = dirname + 'dv'+ '.mat'
       fnameX = dirname + 'xinv'+ '.mat' ! inverse matrix
       fnameK = dirname + 'kv'+ '.mat'
!======================================
       d = spsMat(1,-1,4096,2048,433,433,41,127)
       x(1,1)=d(1,2) 
       x(1,2)=d(1,1)
       if all(abs(d(2:length(d),1))) <tol*ones(fpa_x*fpa_y,1) then
!    % compressed form
          x(2,1)=0
          x(2,2)=-1
       else
          tpv =(d(2:length(d),1)*d(2:length(d),1))
          do j= 2,length(d)
             x(j,1)=d(j,1)/tpv ! % get the scaling
             x(j,2)=d(j,2)
          enddo
       endif
       call fsave(fnameX,x)
    
!       tic; !% start time
       
       if (taskid.eq.MASTER) then
          e=etime(t)
       endif
!============================================
!      Begin main routine
!============================================
       do j=2:n_hcol !vxl_depth %n% j is also pixel id
       write(*,*) 'started column ',j,' | ',e,' sec'
       write(*,*) '--- making d ---'
       d=sps_HxV(fnameX,spsMat(j,-1,4096,2048,433,433,41,127),fnameD,tol)
       write(*,*) '--- making c ---' 
       c = sps_A_MxV(1,(j-1),d,tol);
!     
      if all(abs(c(2:length(c),1))<tol*ones(length(c)-1,1))
        kt=sps_VxH(d,fnameX,(1+TRANSPOSE(d(2:length(d),1))*d(2:length(d),1))); 
!     return compressed vector kt, row vector, compressed columns
        write(*m*) 'not in tolerance'
        call fsave(fnameK,kt)
      else 
!     ========================================
        kt=zeros(length(c),2);
        kt(2:length(c),1)=c(2:length(c),1)/TRANSPOSE((c(2:length(c),1))*c(2:length(c),1)); 
!     % return row vector, compressed rows
        kt(2:length(c),2)=c(2:length(c),2); ! keep the row id's which have become col id's
        kt(1,1)=c(1,2);
        kt(1,2)=c(1,1);
        call fsave(fnameK,kt)
!        msg = 'in tolerance'
      endif
!     use compressed X as well
!     check out nonzero d elements
!     take and save difference of Xprev and new column id in sparse format
!     signed column, row successive
!     convert new kt to signed column, row successive format and append to 
!     fnameX;
    
     call sps_VxVt(fnameD,fnameK,fnameC,tol); 
! output signed column, row successive
     write(*,*) '---- made subtractand ----'
       SELECTCASE(j) !final tolerance setting
          CASE(vxl_depth)
          call sps_sum(fnameX,fnameC,fnameX,-1,ftol) 
          ! return signed column, row-successive format
          CASE DEFAULT
          call sps_sum(fnameX,fnameC,fnameX,-1,tol) 
          ! return signed column, row-successive format
       END SELECT
          CALL X_append(fnameX,kt(2:length(kt),:)); 
          ! snippet off the header for kt
       if (taskid.eq.MASTER) then
          e=etime(t)
       endif
      write(*,*) '========= completed column 'j,' | 'e ' secs ========='
       enddo
       write(*,*) '***** CALCULATED BASIS H PSEUDOINVERSE *****'
      end



!======================================================================
! ROUTINE NAME: part_greville
!
! DESCRIPTION:
!   This routine creates the expanded explicit inverse matrix given
!   H and Hinv in sparse row and sparse col format, respectively.
!   Additinal input is X and Y pixel dimension.  
!
! HISTORY:
!   September 21, 2005: Original Code
!
!   Author:           Latsavongsakda Sethaphong
!   Institution:      Molecular Physiology & Biophysics
!                     Piston Lab
!                     Vanderbilt University
!                     Nashville, TN
!                     latsavongsakda.sethapong@vanderbilt.edu
!----------------------------------------------------------------------
  
       subroutine part_greville()
       end

       subroutine open_data(fname, action)

       end

!======================================================================
! ROUTINE NAME: sps_geninv
!
! DESCRIPTION:
!   This algorithm computes the generalized inverse of any matrix 
!   based upon the Moore Penrose Inverse theorem and the reverse
!   order law
!   Adapted from:
!   Fast Computation of Moore-Penrose Inverse Matrices
!   Pierre Courrieu
!   Neural Information Processing - Letters and Reviews
!   Vol. 8, No.2, August 2005
!
! HISTORY:
!   September 21, 2005: Original Code
!
!   Author:           Latsavongsakda Sethaphong
!   Institution:      Molecular Physiology & Biophysics
!                     Piston Lab
!                     Vanderbilt University
!                     Nashville, TN
!                     latsavongsakda.sethapong@vanderbilt.edu
!----------------------------------------------------------------------

       subroutine sps_geninv(fnameH, fnameO,m,n)
       char*80 fnameH, fnameO
       integer m,n
!      function Y = geninv(G)
!% Returns the Moore-Penrose inverse of the argument
!% Transpose if m < n
       open(unit=10, file='hbasis.dat',status='old')
       i = 1
       do while(.TRUE.)
          read(13,end=45,fmt=*) hbasis(i,1), hbasis(i,2)
          i=i+1
       enddo
45     close(13)

       [m,n]=size(G); transpose=false;
       if m<n
       transpose=true;
       A=G*G';
       n=m;
       else
       A=G'*G;
       end
!  Full rank Cholesky factorization of A
       dA=diag(A); tol= min(dA(dA>0))*1e-9;
       L=zeros(size(A));
       r=0;
       for k=1:n
       r=r+1;
       L(k:n,r)=A(k:n,k)-L(k:n,1:(r-1))*L(k,1:(r-1))';
!% Note: for r=1, the substracted vector is zero
       if L(k,r)>tol
       L(k,r)=sqrt(L(k,r));
       if k<n
       L((k+1):n,r)=L((k+1):n,r)/L(k,r);
       end
       else
       r=r-1;
       end
       end
       L=L(:,1:r);
!% Computation of the generalized inverse of G
       M=inv(L'*L);
       if transpose
          Y=G'*L*M*M*L';
       else
          Y=L*M*M*L'*G';
       endif
       end subroutine
