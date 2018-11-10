function ntG = CGH_sim(src_fname,out_fname,detector_D,separation_d,FPA_sizeX,FPA_sizeY,rad_1,rad_2,rad_3,pxl_cnt, dim_mm_xy);
% give detector_D in millimeters (must be the same as the dimensions on the
% FPA face in pixels per millimeter
% give separation_d in meters
% output of displacement will be in meters (wavelength is in nanometers)
% Create a 7 x 7 array of diffraction orders corresponding to 
% -3,-2,-1,0,1,2,3 in specific radial directions
% 3D Space Mapping and projection onto 2D FPA
% radiance distribution is as follows:
% Itheta = Io * (sin(beta_o)*cos(gamma))^2
% Number of slits = 2
RESX = FPA_sizeX; % number columns
RESY = FPA_sizeY; % number rows
m_orders = 7; % -3,-2,-1,0,1,2,3
srcObj_rgb = imread(src_fname);
bsc = size(srcObj_rgb);
srcObj = rgb2wave(srcObj_rgb); % wave profile
srcObj_Int = rgb2gray(srcObj_rgb); %intensity profile
tgt_sX = bsc(1);
tgt_sY = bsc(2);
%ntD = zeros(FPA_sizeX, FPA_sizeY,1);
% row, column
ntG = zeros(RESY,RESX,1);
%ntG = sparse(RESY, RESX,1);
imwrite(ntG,'image\blank.bmp','bmp');
ntG = imread('image\blank.bmp');
%ntD = imread('blank.bmp');
mapIndex(m_orders^3,2)= double(0.0); % record intensities and remap to a different space
mapInty(m_orders^3,1) = double(0.0);
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
deflection(3,(m_orders-1)) = double(0.0); 
order_m = 1;
omega = 0.001/1200;% meters 600 lines per mm
intys(3,1) = double(0.00000000);
for syp = 1:tgt_sX
    for sxp = 1:tgt_sY
    % establish deflection distance
    % wavelengths 400nm to 800nm
    if srcObj_Int(syp,sxp) > 0
% transform luminesence for three (x',y') -> (x'',y'') -> (x''',y''')
% we get simulate a sequence of three phase gratings rotated in 60 degree
% increments to achieve dispersion in multiple directions and orders
% determine angle sent out first from the normal
% IMPORTANT !!!!!!!!==================== deflection is in millimeters
% millimeters to pixels
mm2pxl_xy = pxl_cnt/dim_mm_xy;
lamda = srcObj(syp,sxp)*10+390; % wave profile map
y = (order_m*lamda*(1e-9)*detector_D/separation_d)*mm2pxl_xy;
mr1 = y;
mr2 = 2*y;
mr3 = 3*y;
% scaled to the intensity at center m=0
theta1 = asin((mr1/mm2pxl_xy)/detector_D); % first order
theta2 = asin((mr2/mm2pxl_xy)/detector_D); % second order
theta3 = asin((mr3/mm2pxl_xy)/detector_D); % third order
%
% G = [cos(t) sin(t); -sin(t) cos(t)]
% deflection from landing on the imaging plane
% first project rotation

beta_1 = pi*omega*(sin(theta1))/(lamda*(1e-9));
beta_2 = pi*omega*(sin(theta2))/(lamda*(1e-9));
beta_3 = pi*omega*(sin(theta3))/(lamda*(1e-9));
gamma_1 = (pi*separation_d/(lamda*(1e-9)))*(sin(theta1));
gamma_2 = (pi*separation_d/(lamda*(1e-9)))*(sin(theta2));
gamma_3 = (pi*separation_d/(lamda*(1e-9)))*(sin(theta3));
% intensity based on projection angle
intys(1) = ((sin(beta_1)/beta_1)*cos(gamma_1))^2;
intys(2) = ((sin(beta_2)/beta_2)*cos(gamma_2))^2;
intys(3) = ((sin(beta_3)/beta_3)*cos(gamma_3))^2;
% rotation angles
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
%mapInty(1)= double(1.0); % start at one
mapInty(1)=double(srcObj_Int(syp,sxp)); % changing to calibration based on image
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
           inM = 1+cq+3*(bp-1);
        if n_Y > 0 && n_Y < RESY + 1 && n_X > 0 && n_X < RESX + 1
           mapIndex(inM,1)=n_Y;
           mapIndex(inM,2)=n_X;
           mapInty(inM)=mapInty(1)*intys(cq);
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
                    inM = (m_orders^ir)+(im-1)*6+iq+3*(ip-1);
                 if n_Y > 0 && n_Y < RESY + 1 && n_X > 0 && n_X < RESX + 1
                    mapIndex(inM,1)=n_Y;
                    mapIndex(inM,2)=n_X;
                    mapInty(inM)=mapInty(inM)+mapInty(im)*intys(cq);
                 end
           end
        end    
    end
end

% mapIndex into the ntD array
cnt = 0;
for iu = 1:(m_orders^3)
        cnt = cnt+1;
        m_Y = int16(mapIndex(iu,1));
        m_X = int16(mapIndex(iu,2));
         m_I = mapInty(iu);
        if m_X < RESX + 1 && m_X > 0 && m_Y < RESY + 1 && m_Y > 0 
           ntG(m_Y,m_X) = double(ntG(m_Y,m_X)) + m_I;% had been reduced
        end
end
% 
        end % for line 49
    end % for line 44
end % for line 43
imwrite(uint8(double(ntG)),out_fname,'bmp');
%catch
%    m_I
%    cnt
    %end