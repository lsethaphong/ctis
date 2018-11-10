function nz_basis = sps_cgh2voxel(Im);
% 14 SEP '04 Latsavongsakda Sethaphong
% RGB space conversion to voxel space
% version 1.0 -- working
% assume we get square pictures
   RES = size(Im);
   nz_pts = 127; % expect no more than this
   col_V(2*nz_pts,1)=0.0;
   high_val = 1.0;
   cnt = 0;
   cnt2 = 0;
   for j=1:RES(1) % row
       for i=1:RES(2) % column
           cnt = cnt+1; % colId
           if Im(j,i) > 0.0
               if Im(j,i) > high_val
               high_val = Im(j,i);
               end
               cnt2 = cnt2 + 1;
               col_V(cnt2,1)= Im(j,i);
               col_V(cnt2+nz_pts,1) = cnt; % keep track of the colId
           end
       end
   end
%   cnt_cgh2voxel = cnt
%   non_zero_cnt = cnt2
   col_V(1:nz_pts,1)= col_V(1:nz_pts,1)/high_val;
   nz_basis=col_V;