% a must be lower triangular and have nxn dimensions
% for use with the geninv.m function
function Lnew=Linv(a)
p = a;
dim = size(a);
      for i=1:dim(1)
         a(i,i)=1./p(i,i);
          for j=(i+1):dim(1)
            sum=0.;
            for k=i:(j-1)
               sum=sum-a(j,k)*a(k,i);
            end
            a(j,i)=sum/p(j,j);
        end 
    end
Lnew=a;
return