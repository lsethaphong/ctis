function to_sgcol(A,fnameA); % signed column, successive row sparse format
%msg = 'signed column, row successive format'
%clear msg;
col=zeros(1,2);
sz=size(A);
cnt = 1;
col(1,1)=sz(1); % record number of rows
col(1,2)=sz(2); % record number of columns
tol = 5.0e-7;
for j=1:sz(1)
%    sign = 1; % beginning new row
    sign = -1;
    for i=1:sz(2)
        if abs(A(j,i)) > tol
           cnt=cnt+1;
           col(cnt,1)=A(j,i);
           col(cnt,2)=i; % record column id
           sign = 1;
        end
    end
%    sign= -1;
    % if we exited without assigning anything then
    if sign == -1
        cnt = cnt +1;
        col(cnt,1)= 0; % special condition for blank row
        col(cnt,2)=-1;
    else
        col(cnt,2)=-col(cnt,2); % negative to signify end of row
    end
end
if fnameA == 'c:/matlab_sv13/work/xcol.mat'
    x=col;
%    clear col;
   save(fnameA,'x','-mat');
else
   save(fnameA,'col','-mat');
end
clear col; % temp show -- comment out later
return;