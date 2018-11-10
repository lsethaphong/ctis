function H_basis = sps_construct_basis(src_size,FPA_sizeX,FPA_sizeY,vxl_depth);
% SPARSE REPRESENTATION of the shift invariant basis matrix
% ntG = H.ntF + N
% 128x128x41 voxels N = 671 744
% g is 671 744 x 1
% f is 171 966 464 x 1
% H is 171 966 464 x 671 744 = 115 517 440 393 216
% step through voxels and record the resulting vector f'
% work on physical dimensions to generate the H transform matrix
%
fname = 'c:/matlab_sv13/work/matdata.dat';
RES = src_size; % 10
nz_pts = 127; % non zero projects per pixel of source
% also depends on the rotation angles involved it is lower than 343 due to
% mapping on the same point
%max = FPA_size; % FPA dimensions
ntG = zeros(FPA_sizeY,FPA_sizeX,1);% create a new tiff FPA target file
ntF = zeros(RES,RES,3);% create false source
imwrite(ntF,'image\blank.bmp','bmp');
ntF = imread('image\blank.bmp','bmp');
con_H(2*nz_pts,vxl_depth)=0.0;
%con_H=sparse(FPA_sizeX*FPA_sizeY,src_size*src_size*vxl_depth); % B & W also voxel space
% read it back in
cnt = 0;
jset = (src_size-mod(src_size,2))/2 + 1;
% take the center point projection
for j=jset:jset % row
    for i=jset:jset %column
        % perform calibration
                cnt = cnt +1; % for the final output
        for idx =1:vxl_depth % for every color for CTIS simulated sensitivity
           % idx
            for channel=1:3
                 % transfer a single wavelength
                ntF(j,i,channel) = wave2rgb(idx,channel);
                % ntG(i,j,1:3) = wave2rgb(idx,1:3);
            end
            imwrite(ntF,'image\ntF_src.bmp','bmp');
            %ntF_gray = rgb2gray(ntF) % convert to gray
            % increasing detector distance to 50 cm
            % was .000014 now .000007
            ntG = CGH_sim('image\ntF_src.bmp','image\ntG_cal.bmp',50,0.000014,FPA_sizeX,FPA_sizeY,60,120,180,512,8);
            % voxel basis vector construction
            % normalize ntG
            c_idx= (cnt-1)*vxl_depth+idx; % voxel number
            % record intensities only.
            %scale_value = 1.0/double(ntF_gray(i,j))
            ntG_dbl = double(ntG); % switch to floating point
            clear ntG;
%            intensity patterns
            con_H(:,c_idx) = sps_cgh2voxel(ntG_dbl); % take only the nonzero portions
            ntF(j,i,:)=uint8(0);
            % need a spectral power relationship
        end
        % reset pixel
    end
end
% achieve single mapping to points on the target FPA
% size(con_H)
%c_idx
save(fname,'con_H','-ASCII');
H_basis = con_H;