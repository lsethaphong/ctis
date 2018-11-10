function vtc = sps_A_MxV(col_start,col_stop,comV,tol); % return compressed row
% 12 OCT 04
% Latsavongsakda Sethaphong
% Performs V-M.comV
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% special call to compressed system matrix
% spsMat(colId,rowId,fpa_x,fpa_y,src_x,src_y,vxl_depth,nz_pts);
% gets returned as signed row
vt = zeros(4096*2048,1); % column summing vector
idx_vt = zeros(1,1); % nonzero indexing vector for row id storage
%tol = 5.0e-7;
idx = 0; % non zero counter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
suband = spsMat(col_stop+1,-1,4096,2048,433,433,41,127);
idx=length(suband)-1;
idx_vt(1:idx)=abs(suband(2:length(suband),2)); % record index
% end of preassignment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for colId=col_start:col_stop
   colV = spsMat(colId,-1,4096,2048,433,433,41,127); % return individual columns, with row id's
   for j=2:length(colV) % basis compressed column vector
      for i=2:length(comV) % compressed column vector
         if colId == abs(comV(i,2))
            % direct indexing of summing vector row elements
            % negative due to subtraction
            vt(abs(colV(j,2)))= vt(abs(colV(j,2))) - colV(j,1)*comV(i,1); % row id
            if abs(vt(abs(colV(j,2)))) > tol
               idx = idx+1;
               idx_vt(idx)=abs(colV(j,2)); % record nonzero row id's
            end
         end
      end
   end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% add V(j,:) back in
for n=2:length(suband)
  vt(abs(suband(n,2))) = vt(abs(suband(n,2))) + suband(n,1); % make suband
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vtc=zeros(1,2);
vtc(1,1)=4096*2048; % number of rows
vtc(1,2)=comV(1,2);
clear colV;
clear comV;
% Special subtraction routine
% sort row id's in ascending order and store in reduced index
idx =1;
rowstore=sort(idx_vt);
clear idx_vt;
prev =0;
for k=1:length(rowstore)
    if prev ~= rowstore(k)
       if abs(vt(rowstore(k))) > tol
          idx = idx+1;
          vtc(idx,1) = vt(rowstore(k)); % record value
          vtc(idx,2) = rowstore(k); % record index
       end
       prev = rowstore(k);
    end
 end
 if idx == 1
   idx=idx+1;
    vtc(idx,2)=-1;
    vtc(idx,1)=0;
 else
    vtc(idx,2)=-vtc(idx,2); % end of row
 end
 return;