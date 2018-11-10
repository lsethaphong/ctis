function v_Im = cgh2voxel(Im);
% 11 AUG '04 Latsavongsakda Sethaphong
% RGB space conversion to voxel space
% version 1.0 -- working
% assume we get square pictures
   RES = size(Im)
   col_V(RES(1)*RES(2),1)=0.0;
   cnt = 0;
   for j=1:RES(1) % row
       for i=1:RES(2) % column
           cnt = cnt+1;
           col_V(cnt,1)=double(Im(j,i))*1.0; 
       end
   end
%   cnt_cgh2voxel = cnt
   v_Im=col_V;
   % save the image as an ascii data file
   save('c:/ctis/comm/g_r.dat','col_V','-ASCII');