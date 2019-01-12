!================================================
! sps_matmul.f
!
! fortran implementation of sparse matrix mult
! 10 NOV '04
!================================================

       program sps_matmul(fnameA,fnameB,
     ~ fnameC,func,b_ki,b_kn,tol) 
       char*80 fnameA, fnameB, fnameC
       real, dimension(,2) :: c_ik
       integer, INTENT(IN) :: func,b_ki,b_kn
       integer nxt_a_ptr, idx
       real, INTENT(IN) ::  tol
       real tmp_sum
 
!================================================
! b_ki, b_kn is only if func is -1
! seems to work 21 SEP 04
! Latsavongsakda Sethaphong
! AxB = C
! aij.bjk=cik
! Using Signed-Row Successive, Sparse Format
! B is sorted into Signed-Column Successive, Sparse Format
! C is output as Signed-Row Successive, Sparse Format
! first entry is ColId = col dimension, val = row dimension
! third entry is the first row's first non zero column entry, with
! load(fnameA,varA);  % signed column, row successive format
! A = x;
! load(fnameB,varB);% signed row, column successive format
! B = h;
! c_ik=zeros(1,2);
! a_ij*b_jk = c_ik
!    A(idx,1) % value % Signed Row Successive Format
!    A(idx,2) % colId 
!    B(idx,1) % value % Signed Column Successive Format
!    B(idx,2) % rowId
!================================================

       nxt_a_ptr = 2
       tmp_sum = 0
       idx = 1
!%c_ik(1,1)=0;
!%c_ik(1,2)=0;
       b_ptr = 2
       col_cnt = 0
       row_cnt = 0
       ik_cnt = 0

!      tol = 5.0e-4
      if (func.eq.1) then ! A is signed column, row successive
                   ! B is signed row, column successive
!    % summing vector
      tmp=zeros(A(1,2),1)
!      idx_vec=zeros(1,1)
!      c_ik=zeros(1,2)
      idx = 1
      sum=0
      do i = 2,length(A)
        tmp(abs(A(i,2))) = A(i,1)
        if A(i,2) < 0
            row_cnt = row_cnt +1      
            do j = 2,length(B)
               sum  = sum + B(j,1)*tmp(abs(B(j,2)))
               if B(j,2) < 0
                  col_cnt = col_cnt + 1
                  if abs(sum) > tol
                     idx=idx+1
                     c_ik(idx,1)= sum
                     c_ik(idx,2)= col_cnt
                  end
                  sum =0
               end
            end
            col_cnt =0
            c_ik(idx,2) = -c_ik(idx,2)! % end of row element
            tmp=zeros(A(1,2),1)
        end
       end
       else !% -1 or something else
!% for two signed column, row successive represented matricies

       if (A(1,2) ~= B(1,1)) ! columns of A must equal
          error_msg =' dimensions not compatible '
          prod = 0
       return
       end
   
       if b_ki < 1
          error_msg = ' invalid start index '
          prod = 0
       return
       end
   
       if b_ki.ge.b_kn
          error_msg=' start index is higher than stop index'
          prod = 0
       return;
       end
! otherwise compute

!  determine length partition of second matrix
       if b_kn.eq.-1
          b_ks = B(1,2)! % take all columns
       else
          b_ks=b_kn
       end
   
       while row_cnt.le.A(1,1)
          row_cnt = 1 + row_cnt; % reset col_cnt
       col_cnt = 0;
       sign = -1;
       a_ptr = 0;
       for b_k=b_ki:b_ks % looping through the columns of B
        col_cnt = col_cnt + 1;
        row_match = 1; % reset row_match
        
       for i = 2:length(B) 
 % "row" loop successively through all the rows of B to
 % for the case that a row has been passed but there has been no match for
 % the current A row's column entry, move on to the next row's column item
 if (A(nxt_a_ptr + a_ptr,2) > 0 && nxt_a_ptr + a_ptr < length(A) && abs(A(nxt_a_ptr + a_ptr,2)) < row_match)
    a_ptr = a_ptr + 1; % must be a unified conditional for the case
 end
          if abs(A(nxt_a_ptr + a_ptr,2)) == row_match && abs(B(i,2)) == b_k % column id
             tmp_sum = tmp_sum + A(nxt_a_ptr + a_ptr,1)*B(i,1);
             if A(nxt_a_ptr + a_ptr,2) > 0 && nxt_a_ptr + a_ptr < length(A) % if not end of row
                a_ptr = a_ptr + 1; % go the next non-zero A entry for this row
            end
          end
      
          if B(i,2) < 0
              % what is the current row
             row_match = row_match + 1; % moved past a row
          end
         
      end  % end of "row" loop   
% sum at the end of looping through row for given column k of b
         if abs(tmp_sum) > 0
            idx=idx+1; % increment to next result element
            c_ik(idx,1) = tmp_sum;
            c_ik(idx,2) = b_k; % col id
            tmp_sum = 0; % reset tmp_sum
            sign = 0;
         end
         if b_k < b_ks
            a_ptr = 0; % reset to beginning of current row
         end
         %     end
      
   end % end A-B row loop
       if sign < 0
         idx=idx+1;
         c_ik(idx,1)=0; % blank entry
         c_ik(idx,2) = -1;
       else
         c_ik(idx,2)=-c_ik(idx,2);
       end
       
      if nxt_a_ptr + a_ptr < length(A)
         nxt_a_ptr = nxt_a_ptr + a_ptr +1;% move
         a_ptr = 0;
      end
      
   end % end A row loop loop
end
% record results
c_ik(1,1)=row_cnt;
c_ik(1,2)=B(1,2);
prod=c_ik
save(fnameC,'c_ik','-ASCII');
clear c_ik;
return;
