function pwr_dist_V = map2pwrdist(ntF,img_typ,src_Int);
   RES = length(ntF);
   col_V(RES*RES,1)=0.0; % returns a column vector
   tmp_ntF = ntF;
   if img_typ == 1
   tmp_ntF=rgb2gray(ntF);
   end
   % intensity mapping only -- no voxels
   % monochromatic with intensity values
   for i=1:RES
       for j=1:RES
           col_V(i*j)= double(tmp_ntF(i,j))/src_Int;
       end    
   end
   pwr_dist_V = col_V;