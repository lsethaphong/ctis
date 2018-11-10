function newsim;
Hbasis=sps_construct_basis(433,4096,2048,41);
clear all;
fim=imread('c:/matlab_sv13/work/image/tst433_6.bmp');
f_vxl = map2voxel(fim,41);% must be 433by433
clear fim;
sps_CGH_shift(f_vxl,4096,2048,433,433,41,127,'c:/matlab_sv13/work/image/spsCGH.bmp');