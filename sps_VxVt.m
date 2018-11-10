function sps_VxVt(fnameA,fnameB,fnameC,tol);
% seems to work 21 SEP 04
% Latsavongsakda Sethaphong
% AxB = C
% aij.bjk=cik
% A is taken as Signed-Row, Column Successive, Sparse Format
% B is taken as Signed-Column, Row Successive, Sparse Format
% C is output as Signed-Column, Row-Successive, Sparse Format
% first entry is ColId = col dimension, val = row dimension
% third entry is the first row's first non zero column entry, with
% abs(ColId) being the col Id
%fid = fopen(fnameA);
%%%%%[a] = fscanf(fid,'%e',[vxl_depth nz_pts*2]);
%[a] = fscanf(fid,'%g',[2 inf]); % two rows
%A = a';  % signed row, col successive
%fclose(fid);
load(fnameA,'d');
%A = d;
%clear d;
%size(A)
%clear a;
%fid2 = fopen(fnameB);
%[b] = fscanf(fid2,'%g',[2 inf]);
%fclose(fid2);
%B = b';% signed column, row successive
load(fnameB,'kt');
%B= kt;
%clear kt;
%size(B)
%clear b;
%nz_pts_b = size(B); sort into rows
%nz_pts_a = size(A);
cik=zeros(1,2);
cik(1,1)=d(1,1);
cik(1,2)=kt(1,2);
% a_ij*b_jk = c_ik
%for idx = 1:nz_pts % cycling through
%    A(idx,1) % value % Signed Row Successive Format
%    A(idx,2) % colId 
%    B(idx,1) % value % Signed Column Successive Format
%    B(idx,2) % rowId
    %loop through rows of A
%ik_cnt = 1;
% creates a sparse matrix representation of X output in signed column row
% successive format with rowskip for large zero rows.
ik_cnt = 1;
%tol=5.0e-7;
for i=2:length(d)
    for j=2:length(kt)
        tmp = d(i,1)*kt(j,1);
        if abs(tmp) > tol % it must be greater than tolerance
          ik_cnt = ik_cnt + 1;
          cik(ik_cnt,1) = tmp;% new value
          cik(ik_cnt,2) = abs(kt(j,2));% column id
        end
    end
    cik(ik_cnt,2) = -cik(ik_cnt,2); % end of the row, hence negative column id
end
save(fnameC,'cik','-mat'); % output as signed column, row successive format
clear cik;
return;
