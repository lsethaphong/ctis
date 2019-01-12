       function sps_HxV(fnameX,din,fnameD,tol); %
!% Latsavongsakda Sethaphong 
!% this function multiplies the first few columns of 
!% the basis set by the growing V as part of the inversion
!% algorithm
!% returns a signed column, row successive format
!% H is signed column row successive.
!% H is condensed to signed col, row successive format
!% 
! this program can be expanded to generalized form into 
! hpinv HxV
load(fnameX,'x');
if x(1,2) ~= din(1,1)
    msg ='------- column x not equal to rows of din ---------'
    return;
else
   d=zeros(1,2);
   ik = 1;
   rowid = 0;
   d(1,1)=x(1,1);
   d(1,2)=1;
   xvec = zeros(x(1,2),1); !% number of rows
   summand =0;    
   for j =2:length(x)
      xvec(abs(x(j,2))) = x(j,1); !% take the value
      if x(j,2) < 0 !% reached end of row element
         rowid = rowid + 1;
              for i=2:length(din) !% loop through all elements of d
              summand = summand + din(i,1)*xvec(abs(din(i,2)));
           end
           if abs(summand) > tol
              ik=ik+1;
              d(ik,1)=summand;
              d(ik,2)=rowid;
           end
           xvec=zeros(x(1,2),1); !% clear xvec for next element
           summand = 0;
       end
       end
       d(ik,2)=-d(ik,2); % end of the column
       save(fnameD,'d','-mat');
       return;
       end


