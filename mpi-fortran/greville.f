!==============================================
!    
!     greville.f
!
!     translated by:
!     Latsavongsakda Sethaphong, 7 SEP 04
!     VUMC, Piston Lab
!
!===============================================

        subroutine GREVILLE(A, mi, ni, Hp)
        integer :: i,ii,iii,j,k,mi,ni
        dimension d(mi,1), c(mi,1), tmp_d(1,1),val_d(1,1), v_one(1,1)
        dimension X(ni,mi),dt(1,ni), bt(1,mi),tmp_X(ni,mi)
        dimension A(mi,ni)
        dimension Hp(ni,mi) 
        real max_d,min_d,max_c,min_c,tol,val_d,v_one
        tol = 5.0e-7 
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
               dt=TRANSPOSE(d)
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
