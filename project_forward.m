function trIm = project_forward(Fim,Hmatrix)
%
%
% convert source image to voxel space
v_Fim = 
pfIm = Hmatrix*v_Fim;
% convert pfIm from voxel space to RGB image space
