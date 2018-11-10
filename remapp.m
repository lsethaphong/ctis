function remapp(srcfile,outfile,vxl_depth)
iMap2 = imread(srcfile);
vMapp2 = map2voxel(iMap2,vxl_depth);
voxel2rgb(outfile,vMapp2,vxl_depth);