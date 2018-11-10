function recon4voxel;
fname = 'c:\matlab_sv13\work\image\f_r';
fid = fopen('c:\ctis\comm\f_f.dat');
[a] = fscanf(fid,'%e',[inf]); % It has two rows now.
fclose(fid);
%voxel2rgb(fname,a,41);
full2channels(fname,a,280);
clear a;
%imshow(imread(fname));
return