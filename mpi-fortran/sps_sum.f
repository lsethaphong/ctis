
       subroutine sps_sum(fnameA,fnameB,fnameS,opera,tol);
! 12 NOV 05
! Latsavongsakda Sethaphong
! A and B should have the same dimensions
! then, A +/- B => C requires only row and column align on 
! successive rows with signed column elements
! format is signed column, row successive
       file 
 load(fnameA,'x');
A = x;
clear x;
load(fnameB,'cik'); ! VxVt output
B = cik; ! multiplicative output
clear cik; 
! -------------------------------------------------------
! commit to operations
! -------------------------------------------------------
n_x = zeros(1,A(1,2)); ! row vector
b_ptr = 2;
x = zeros(1,2); ! start off with zero rows
x(1,1)=0; ! number of initial rows
x(1,2)=A(1,2); ! number of expected columns
save(fnameS,'x'); %replace x
%tol = 5.0e-7;
if B(1,1) ~= A(1,1) && B(1,2) ~= A(1,2)
    error_msg = 'dimensions do not match'
    clear all;
    return;
elseif operator == -1

    ! make col index vector
     idx_vc=zeros(1,1);
     idx_cnt = 0; ! keeps all non zero col id's with idx_vc
    for j = 2:length(A)
        n_x(abs(A(j,2)))=A(j,1); ! index of new difference vector
        idx_cnt = idx_cnt+1;
        idx_vc(idx_cnt)=abs(A(j,2)); ! take the colid's of A
        if A(j,2) < 0
            ! create row vector of B
            while B(b_ptr,2) > 0
               n_x(abs(B(b_ptr,2))) = n_x(abs(B(b_ptr,2))) - B(b_ptr,1);
               idx_cnt = idx_cnt+1;
               idx_vc(idx_cnt)=abs(B(b_ptr,2)); ! take the colid's of B
               b_ptr = b_ptr + 1;
            end
            n_x(abs(B(b_ptr,2))) = n_x(abs(B(b_ptr,2))) - B(b_ptr,1);
            idx_cnt = idx_cnt+1;
            idx_vc(idx_cnt)=abs(B(b_ptr,2)); ! take the colid's of B
            b_ptr = b_ptr + 1;  ! when greater than length(B), reached end
            idx_vc = sort(idx_vc); ! sort ascending
            ! eval and store difference above tol
            prev=0;
            com_cnt=0;
            for q=1:length(idx_vc)
                if prev ~= idx_vc(q)
                   prev=idx_vc(q);
                   if abs(n_x(idx_vc(q))) > tol
                       com_cnt =com_cnt+1;
                       comV(com_cnt,1) = n_x(idx_vc(q)); ! value
                       comV(com_cnt,2) = idx_vc(q); ! col index
                   end
                end
            end
            comV(com_cnt,2)=-comV(com_cnt,2); ! endof column
            ! append row by row
            X_append(fnameS,comV);
            comV=zeros(1,2); ! clear comV
            n_x = zeros(1,A(1,2));
            idx_vc = zeros(1,1); ! recompress index vector
            idx_cnt = 0; ! keeps the idx_vc from growing too long
        end
        if b_ptr > length(B) ! reached the end of B's then
            return; ! we're done
        end
    end
    clear n_x;
elseif operator == 1
      % create row vector of A
     % make col index vector
     idx_vc=zeros(1,1);
     idx_cnt = 0;
    for j = 2:length(A)
        n_x(abs(A(j,2)))=A(j,1); % index of new difference vector
        idx_cnt = idx_cnt+1;
        idx_vc(idx_cnt)=abs(A(j,2)); % take the colid's of A
        if A(j,2) < 0
            % create row vector of B
            while B(b_ptr,2) > 0
               n_x(abs(B(b_ptr,2))) =n_x(abs(B(b_ptr,2)))- B(b_ptr,1);
               idx_cnt = idx_cnt+1;
               idx_vc(idx_cnt)=abs(B(b_ptr,2)); % take the colid's of B
               b_ptr = b_ptr + 1;
            end
            n_x(abs(B(b_ptr,2))) =n_x(abs(B(b_ptr,2)))- B(b_ptr,1); % is negative
            idx_cnt = idx_cnt+1;
            idx_vc(idx_cnt)=abs(B(b_ptr,2)); % take the colid's of B
            b_ptr = b_ptr + 1;  
            idx_vc = sort(idx_vc);
            prev=0;
            com_cnt=0;
            for q=1:length(idx_vc)
                if prev ~= idx_vc(q)
                   prev=idx_vc(q);
                   if abs(n_x(idx_vc(q))) > tol
                       com_cnt =com_cnt+1;
                       comV(com_cnt,1) = n_x(idx_vc(q)); % value
                       comV(com_cnt,2) = idx_vc(q); % col index
                   end
                end
            end
            comV(com_cnt,2)=-comV(com_cnt,2); % end of column
            X_append(fnameS,comV);
            %append row by row
            comV=zeros(1,2); % clear comV
            n_x = zeros(1,A(1,2));            
            idx_vc = zeros(1,1); % recompress index vector
            idx_cnt = 0;
        end
        if b_ptr > length(B) % reached the end of B's then
            return;
        end
    end
    clear n_x;
%    clear bkt;
else
    error_msg = 'unknown operator';
    return;
end
%msg = 'done'
%save(fnameS,'x','-mat');
return;

       end
