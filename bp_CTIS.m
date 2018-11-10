function bp_Im = bp_CTIS(src_fname,out_fname,bp_fname);
% simple test case for BW CTIS to spectral mapping
%construct_basis(src_size,FPA_size,vxl_depth)
Ho=construct_basis(4,64,20);
Hp=pinv(Ho);
f_hat=imread(src_fname);
f_hat
g_hat=CGH_sim(src_fname,out_fname,3.3,0.000014,64,64,60,120,180,512,8);
%g_o = double(g_hat);
g_o = cgh2voxel(g_hat);
%size(g_o)
f_o = Hp*g_o;
%size(f_o)
% reconstruct f_hat to rgb image
%voxel2rgb(fname,col_V,vxl_depth)
%mart(ntG_prior, ntF_prior, Hmat, psinv_Hmat,cycle_stop);
recon_Im=mart(g_o,f_o,Ho,Hp,8,1);
bp_Im=voxel2rgb(bp_fname,recon_Im,20)
