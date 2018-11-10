function ntG = sps_CGH_shift(f_vxl,fpa_x,fpa_y,src_x,src_y,vxl_depth,nz_pts,out_fname);
% f_vxl is a voxel representation of the source image
%try 
ntG(fpa_x*fpa_y)=0.0;
fid = fopen('Ho_sps.dat'); % for the maximum basis set with the modeled code
[a] = fscanf(fid,'%e',[vxl_depth nz_pts*2]); % It has two rows now.
fclose(fid);
basisH = a'; % 254 rows, 41 columns = voxel depth
clear a;
sz_H = size(basisH);
vxl_depth = sz_H(2);
% assume that the source distance is odd
ref_x = (src_x-1)/2 + 1;
ref_y = (src_y-1)/2 + 1;
% the colId gives us the voxel of interest
for colId = 1:length(f_vxl)
   % looping through all of the voxels in the vector source 
   if (f_vxl(colId) > 0) % not wasting time on zero elements
       
vxl_basis_id = mod(colId,vxl_depth);

if vxl_basis_id == 0
    vxl_basis_id = vxl_depth;
end
% determine source displacement
pixel_x_id = mod((colId-mod(colId,vxl_depth)),vxl_depth*src_x)/vxl_depth + 1; % get the pixel 
pixel_y_id = (colId-mod(colId,vxl_depth*src_x))/(vxl_depth*src_x) + 1; % take the whole number
% id's of zero means it's at the bottom right corner or at the last column
% displacement vectors
shift_x = pixel_x_id-ref_x;
shift_y = (pixel_y_id-ref_y)*fpa_x;
%shift_f = shift_x + shift_y*fpa_x; % shifting along the pixels of the FPA which are row elements of the

for ptr = 1:sz_H(1)/2 % goes through the entire listing of elements
    pixel_id=basisH(ptr+sz_H(1)/2,vxl_basis_id);
    x_rem = rem(pixel_id,fpa_x);
    pixel_id = pixel_id + shift_x + shift_y;
%     can't move more than allowed y where pixel_id is [0,4k.2k]     
   if pixel_id <= length(ntG) && pixel_id > 0 && (shift_x + x_rem > 0) && (shift_x + x_rem < fpa_x)
      ntG(pixel_id) = ntG(pixel_id)+ f_vxl(colId)*basisH(ptr,vxl_basis_id);
   end
end
% transform matrix H
% set shift conditions for the rowId's
%t_ptr = 0;
%for ptr = 1:sz_H(1)/2
%   newrow = basisH(ptr+sz_H(1)/2,vxl_basis_id) + shift_f;
%   if rowId == newrow
%       t_ptr = ptr;
%       ptr = sz_H(1)/2; % stop the for loop
%   end
%end
%   if t_ptr > 0
%       valM=basisH(t_ptr,vxl_basis_id);
%   else
%       valM = 0;
%   end
   end
end

%catch
%    vxl_basis_id
%    colId
%    pixel_id
%    x_rem
%    shift_x
%    shift_y
%    ptr
%end

voxel2cgh(ntG,fpa_x,fpa_y);
%imwrite(uint8(double(ntG)),out_fname,'bmp');
