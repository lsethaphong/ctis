function create_sgrow_H(fnameA); % signed column, successive row sparse format
%spsMat(colId,rowId,fpa_x,fpa_y,src_x,src_y,vxl_depth,nz_pts)
msg = 'signed column, row successive format'
clear msg;
col=zeros(1,2);
cnt = 1;
col(1,1)=4096*2048; % record number of rows
col(1,2)=433*433*41; % record number of columns
for j=1:col(1,1)
%    sign = 1; % beginning new row
    sign = -1;
    for i=1:col(1,2) % number of columns
        if spsMat(j,i,4096,2048,433,433,41,127) ~= 0
           cnt=cnt+1;
           col(cnt,1)=spsMat(j,i,4096,2048,433,433,41,127);
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
save(fnameA,'col','-ASCII');
%col % temp show -- comment out later
return;