function full2channels(fname,col_V,vxl_depth);
% 7 MAR '05 Latsavongsakda Sethaphong
% Voxel space to RGB reconstruction algorithm
% version 1.0 -- working
% assume we get square images to reconstruct
%RES = (length(col_V)/vxl_depth)^0.5 % gives the size
objI=zeros(80,89,3); % not square
imwrite(objI,'image\blank.bmp','bmp');
  tmp_h = 0.0039;%double(0.0);
  tmp = double(0.0);
  for channel = 1:vxl_depth
      myname = ['image\f_r' int2str(channel) '.bmp'] % show the name being worked on
      objI=imread('image\blank.bmp');
     for j = 1:89 % row 'y'
        for i = 1:80 % column 'x'
              % needed the pixel id
              tmp = col_V(((j-1)*RES + i-1)*vxl_depth  + channel ); % record position
%              for ck =1:3
%              objI(j,i,:) = uint8(double(tmp*wave2rgb(channel,ck))); % comes back as 1 %record strongest intensity
              %              end
              objI(j,i,1)=uint8(tmp)%*255.0); % when not scaled
          end % end channels
      end
   imwrite(objI,myname,'bmp');
   end
   %-------------------------
return;