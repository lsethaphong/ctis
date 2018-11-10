function rctIm = recoverI(r_file,fname,FPA_sizeX,FPA_sizeY,vxl_depth,obj_size);
% test the reconstruction algorithms and subroutines

  gF = imread(r_file);
  RES = length(gF);
  Ho = construct_basis(obj_size,FPA_sizeX,FPA_sizeY,vxl_depth);
%  Hp = pinv(Ho);
  Hp = greville(Ho);
  % voxel space column vectors
  col_V(RES*RES*vxl_depth,1)=0;
  vxl_rgb(4*4*vxl_depth,1)=0;
  % converting to voxel space
  col_V = map2voxel(gF,vxl_depth); % check out mappings
  % using Moore-Penrose pseudoinverse method
  vxl_rgb = Hp*(Ho*col_V); % conversion of image to voxel space projection and back again
  %vxl_rgb % still in voxel space projections
  objRES = length(vxl_rgb)/(vxl_depth*obj_size);
  ntG_sim = CGH_sim(r_file,'image\tgt11x11_nH_ntG.bmp',10,0.000014,FPA_sizeX,FPA_sizeY,60,120,180,512,8);
  ntG_dbl = (double(ntG_sim));
  clear ntG_sim;
  ntG_v =cgh2voxel(ntG_dbl);
  clear ntG_dbl;
  ntF_v = Hp*ntG_v;
%  mart(ntG_prior, ntF_prior, Ho, Hp,cyclestop);
  ntF_v = mart(ntG_v, ntF_v, Ho, Hp, 7); % simple reconstruction
  clear Hp;
  tmp_V=Ho*col_V;
  voxel2cgh(tmp_V,FPA_sizeX,FPA_sizeY); % creates a grayscale image of projection
  clear Ho;
  clear col_v;
  imwrite(zeros(objRES,objRES,3), fname, 'bmp');
%  objI = imread(fname); 
  % convert voxels to RGB
  message = 'inversion from matrix projection'
  voxel2rgb(fname,vxl_rgb,vxl_depth);
  message = 'direct inversion from computed projection'
  voxel2rgb('image\tgt_ntFv.bmp',ntF_v,vxl_depth);
  %  voxel2rgb(fname,col_V,vxl_depth);

rctIm= imread(fname);
  imshow(rctIM);