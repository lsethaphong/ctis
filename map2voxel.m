function v_Im = map2voxel(Im,vxl_depth);
% 11 AUG '04 Latsavongsakda Sethaphong
% RGB space conversion to voxel space
% version 1.0 -- working
% assume we get square pictures
   RES = size(Im);
   col_V(RES(1)*RES(2)*vxl_depth,1)=0;
   mapF = rgb2wave(Im); % gives the voxel mapping per grid 
   cnt = 0;
   for j=1:RES % row
       for i=1:RES % column
           cnt = cnt+1;
           if mapF(j,i) > 0 & mapF(j,i) < (vxl_depth+1)
           tmp= double(Im(j,i,1))+double(Im(j,i,2))+double(Im(j,i,3));
           tmp= tmp/(double(wave2rgb(mapF(j,i),1))+double(wave2rgb(mapF(j,i),2))+double(wave2rgb(mapF(j,i),3)));
           col_V((cnt-1)*vxl_depth + mapF(j,i))= tmp;
           end    
       end
   end
   %cnt
   v_Im = col_V;