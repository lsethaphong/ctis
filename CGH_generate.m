function ntG = CGH_generate(src_fname,out_fname,detector_D,separation_d,FPA_sizeX,FPA_sizeY,rad_1,rad_2,rad_3);
% give detector_D in millimeters
% give separation_d in meters
% output of displacement will be in meters (wavelength is in nanometers)
% Create a 7 x 7 array of diffraction orders corresponding to 
% -3,-2,-1,0,1,2,3 in specific radial directions
% 3D Space Mapping and projection onto 2D FPA
%
try
RESX = FPA_sizeX;
RESY = FPA_sizeY;
m_orders = 7; % -3,-2,-1,0,1,2,3
%
%ntF = zeros(FPA_size,FPA_size); % maps to a 2D intensity profile
%srcObj = zeros(source_size,source_size,1);
%imwrite(srcObj,'bsrc.bmp','bmp');
srcObj_rgb = imread(src_fname);
bsc = size(srcObj_rgb);
srcObj = rgb2wave(srcObj_rgb); % wave profile
srcObj_Int = rgb2gray(srcObj_rgb); %intensity profile
tgt_sX = bsc(1);
tgt_sY = bsc(2);
ntD = zeros(FPA_sizeX, FPA_sizeY,1);
ntG = zeros(FPA_sizeX, FPA_sizeY,1);
imwrite(ntD,'blank.bmp','bmp');
ntG = imread('blank.bmp');
ntD = imread('blank.bmp');
mapIndex(1:(m_orders^3),1:3)=double(0); % record intensities and remap to a different space
% y ~ m_order*lamda*screen_D/separation_d
% calculate for as many orders as possible to get full coverage of the FPA
% 29 incremental 10nm from 430nm to 710nm
% 41 incremental from 400nm to 800nm
% grating spacings d_1_2, d_2_3, d_3_FPA
% let d_1_2 = d_2_3 = 1 cm
% let d_3_FPA = 10 cm
% let the gratings be a 5cmx5cm 
% let the FPA dimensions be 30cmx30cm 
origin_X = FPA_sizeX/2;
origin_Y = FPA_sizeY/2;
origin_Xs = tgt_sX/2;
origin_Ys = tgt_sY/2;
deflection = zeros(3,m_orders-1); 
order_m = 1;
cnt = 0;
cnt1 = 0;
cnt2 = 0;
cnt3 = 0;
for sxp = 1:tgt_sX
    for syp = 1:tgt_sY
%        for lamda =1:vxl_depth %41
%            srcObj(sxp,syp) = 255; % intensity only
    % establish deflection distance
    % wavelengths 400nm to 800nm
    if srcObj_Int(sxp,syp) > 0
    lamda = srcObj(sxp,syp); % wave profile map
    y = order_m*(390+10*lamda)*(1e-9)*detector_D/separation_d;
% transform luminesence for three (x',y') -> (x'',y'') -> (x''',y''')
% we get simulate a sequence of three phase gratings rotated in 60 degree
% increments to achieve dispersion in multiple directions and orders
% determine angle sent out first from the normal
mr1 = y;
mr2 = 2*y;
mr3 = 3*y;
theta1 = asin(mr1/detector_D); % first order
theta2 = asin(mr2/detector_D); % second order
theta3 = asin(mr3/detector_D); % third order
% G = [cos(t) sin(t); -sin(t) cos(t)]
% deflection from landing on the imaging plane
% first project rotation
phi_1 = (rad_1/180)*pi;
phi_2 = (rad_2/180)*pi;
phi_3 = (rad_3/180)*pi;

deflection(1,1)=mr1*sin(phi_1);
deflection(1,3)=mr2*sin(phi_1);
deflection(1,5)=mr3*sin(phi_1);
deflection(1,2)=mr1*cos(phi_1);
deflection(1,4)=mr2*cos(phi_1);
deflection(1,6)=mr3*cos(phi_1);
% second projection rotation
deflection(2,1)=mr1*sin(phi_2);
deflection(2,3)=mr2*sin(phi_2);
deflection(2,5)=mr3*sin(phi_2);
deflection(2,2)=mr1*cos(phi_2);
deflection(2,4)=mr2*cos(phi_2);
deflection(2,6)=mr3*cos(phi_2);
% third projection rotation
deflection(3,1)=mr1*sin(phi_3);
deflection(3,3)=mr2*sin(phi_3);
deflection(3,5)=mr3*sin(phi_3);
deflection(3,2)=mr1*cos(phi_3);
deflection(3,4)=mr2*cos(phi_3);
deflection(3,6)=mr3*cos(phi_3);
% project onto 2D FPA space some distance away.
%end
n_Y =origin_Y + syp-origin_Ys;
n_X =origin_X + sxp-origin_Xs;
mapIndex(1,1)=n_Y;
mapIndex(1,2)=n_X;
mapIndex(1,3)=double(srcObj_Int(sxp,syp));
for bp = 1:2
    for cq = 1:((m_orders-1)/2)
        i_y = 1+(cq-1)*2;
        i_x = 2+(cq-1)*2;
        if bp == 1 
           n_Y = mapIndex(1,1)+deflection(1,i_y);
           n_X = mapIndex(1,2)+deflection(1,i_x);
        else
           n_Y = mapIndex(1,1)-deflection(1,i_y);
           n_X = mapIndex(1,2)-deflection(1,i_x);
        end
        if n_Y > 0 & n_Y < RESY + 1 & n_X > 0 & n_X < RESX + 1
           mapIndex(1+cq+3*(bp-1),1)=n_Y;
           mapIndex(1+cq+3*(bp-1),2)=n_X;
           mapIndex(1+cq+3*(bp-1),3)=mapIndex(1,3);
       end
    end
end
    % building the full mapping (each point maps onto itself and six other
    % points
for ir = 1:2
    rt = 1+ir;
    im_end = m_orders^ir;
    for  im = 1:im_end %gets us to the 343 final points of projection
        for ip = 1:2
             for iq = 1:((m_orders-1)/2)
                 i_y = 1+(iq-1)*2;
                 i_x = 2+(iq-1)*2;
                 if ip == 1 
                   n_Y = mapIndex(im,1)+deflection(rt,i_y);
                   n_X = mapIndex(im,2)+deflection(rt,i_x);
                 else
                   n_Y = mapIndex(im,1)-deflection(rt,i_y);
                   n_X = mapIndex(im,2)-deflection(rt,i_x);
                 end
                 if n_Y > 0 & n_Y < RESY + 1 & n_X > 0 & n_X < RESX + 1
                    mapIndex((m_orders^ir)+(im-1)*6+iq+3*(ip-1),1)=n_Y;
                    mapIndex((m_orders^ir)+(im-1)*6+iq+3*(ip-1),2)=n_X;
                    mapIndex((m_orders^ir)+(im-1)*6+iq+3*(ip-1),3)=mapIndex(im,3);
                 end
           end
        end    
    end
end

% mapIndex into the ntD array

for iu = 1:(m_orders^3)
    if mapIndex(iu,3) > 0
        cnt = cnt+1;
        m_Y = int16(mapIndex(iu,1));
        m_X = int16(mapIndex(iu,2));
        m_I = int16(mapIndex(iu,3));
        if m_X < RESX + 1 & m_X > 0 & m_Y < RESY + 1 & m_Y > 0 
           ntG(m_X,m_Y) = uint8(double(ntG(m_X,m_Y)) + double(m_I));
        end
    end
end
% 
        end % for line 49
    end % for line 44
end % for line 43
%cnt
imwrite(ntG,out_fname,'bmp');
catch
    m_I
    cnt
end