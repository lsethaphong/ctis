function ntG = CGH_rayleigh(detector_D,separation_d,FPA_sizeY,FPA_sizeX,tgt_sX,tgt_sY,phi_1,phi_2,phi_3,vxl_depth)
% Create a 7 x 7 array of diffraction orders corresponding to 
% -3,-2,-1,0,1,2,3 in specific radial directions
% 3D Space Mapping and projection onto 2D FPA
%
source_size = 128;
RESX = FPA_sizeX;
RESY = FPA_sizeY;
%
%ntF = zeros(FPA_size,FPA_size); % maps to a 2D intensity profile
srcObj = zeros(source_size,source_size,1);
imwrite(srcObj,'bsrc.bmp','bmp');
srcObj = imread('bsrc.bmp');
ntD = zeros(FPA_sizeX, FPA_sizeY,1);
imwrite(ntD,'blank.bmp','bmp');
ntD = imread('blank.bmp');
ntD2 = imread('blank.bmp');
ntD3 = imread('blank.bmp');
% y ~ m_order*lamda*screen_D/separation_d
% calculate for as many orders as possible to get full coverage of the FPA
% 29 incremental 10nm from 430nm to 710nm
% 41 incremental from 400nm to 800nm
% grating spacings d_1_2, d_2_3, d_3_FPA
% let d_1_2 = d_2_3 = 1 cm
% let d_3_FPA = 10 cm
% let the gratings be a 5cmx5cm 
% let the FPA dimensions be 30cmx30cm 
order_m = 1; % initial
middleX = FPA_sizeX/2;
middleY = FPA_sizeY/2
origin_X = middleX;
origin_Y = middleY;
origin_Xs = tgt_sX/2;
origin_Ys = tgt_sY/2;
for sxp = 1:tgt_sX
    for syp = 1:tgt_sY
        for lamda =1:vxl_depth %41
            srcObj(sxp,syp) = 255; % intensity only
    % establish deflection distance
    % wavelengths 400nm to 800nm
    y = order_m*(390+10*lamda)*(1e-9)*detector_D/separation_d
% transform luminesence for three (x',y') -> (x'',y'') -> (x''',y''')
% we get simulate a sequence of three phase gratings rotated in 60 degree
% increments to achieve dispersion in multiple directions and orders
% determine angle sent out first from the normal
mr1 = y;
mr2 = 2*y;
mr3 = 3*y;
theta1 = asin(mr1/detector_D) % first order
theta2 = asin(mr2/detector_D) % second order
theta3 = asin(mr3/detector_D) % third order
% G = [cos(t) sin(t); -sin(t) cos(t)]
% deflection from landing on the imaging plane
% first project rotation = 60 degrees
m1_yp = mr1*sin(phi_1);
m2_yp = mr2*sin(phi_1);
m3_yp = mr3*sin(phi_1);
m1_xp = mr1*cos(phi_1);
m2_xp = mr2*cos(phi_1);
m3_xp = mr3*cos(phi_1);
% generate detection profile
% second projection rotation = 120 degrees of signals detected
m1_ypp = mr1*sin(phi_2);
m2_ypp = mr2*sin(phi_2);
m3_ypp = mr3*sin(phi_2);
m1_xpp = mr1*cos(phi_2);
m2_xpp = mr2*cos(phi_2);
m3_xpp = mr3*cos(phi_2);
% generate detection profile
% third projection rotation = 180 degrees of signals detected
m1_yppp = mr1*sin(phi_3);
m2_yppp = mr2*sin(phi_3);
m3_yppp = mr3*sin(phi_3);
m1_xppp = mr1*cos(phi_3);
m2_xppp = mr2*cos(phi_3);
m3_xppp = mr3*cos(phi_3);
% project onto 2D FPA space some distance away.
end
% 60 degree rotation
for vxp = 1:tgt_sX
    for vyp = 1:tgt_sY
        if(srcObj(vxp,vyp) > 0)
            % redraw coordinates for vxp,vyp
            vxpp = origin_X + vxp-origin_Xs;
            vypp = origin_Y + vyp-origin_Ys;
           ntD(vxpp,vypp)=srcObj(vxp,vyp); % 0
           if ((vxpp+m1_xp) <= RESX) & ((vypp+m1_yp) <= RESY) & (vxpp+m1_xp)>= 1 & (vypp+m1_yp)>= 1
           ntD(uint16(vxpp+m1_xp),uint16(vypp+m1_yp))=srcObj(vxp,vyp); % 1
           end
           if ((vxpp+m2_xp) <= RESX) & ((vypp+m2_yp) <= RESY) & (vxpp+m2_xp)>= 1 & (vypp+m2_yp)>= 1
           ntD(uint16(vxpp+m2_xp),uint16(vypp+m2_yp))=srcObj(vxp,vyp); % 2
           end
           if ((vxpp+m3_xp) <= RESX) & ((vypp+m3_yp) <= RESY) & (vxpp+m3_xp)>= 1 & (vypp+m3_yp)>= 1
           ntD(uint16(vxpp+m3_xp),uint16(vypp+m3_yp))=srcObj(vxp,vyp); % 3
           end
           if ((vxpp-m1_xp) <= RESX) & ((vypp-m1_yp) <= RESY) & (vxpp-m1_xp)>= 1 & (vypp-m1_yp)>= 1
           ntD(uint16(vxpp-m1_xp),uint16(vypp-m1_yp))=srcObj(vxp,vyp); % -1
           end
           if ((vxpp-m2_xp) <= RESX) & ((vypp-m2_yp) <= RESY) & (vxpp-m2_xp)>= 1 & (vypp-m2_yp)>= 1
           ntD(uint16(vxpp-m2_xp),uint16(vypp-m2_yp))=srcObj(vxp,vyp); % -2
           end
           if ((vxpp-m3_xp) <= RESX) & ((vypp-m3_yp) <= RESY) & (vxpp-m3_xp)>= 1 & (vypp-m3_yp)>= 1
           ntD(uint16(vxpp-m3_xp),uint16(vypp-m3_yp))=srcObj(vxp,vyp); % -3       
           end
       end
    end
end
 % blank source
 imwrite(ntD,'CGH3_ntD1.bmp','bmp');
 imwrite(srcObj,'CGH3_ntD0.bmp','bmp');
 srcObj(sxp,syp) = 0;
% 120 degree rotation
for vxp = 1:RESX
    for vyp = 1:RESY
        if(ntD(vxp,vyp) >= 1)
           ntD2(vxp,vyp)=ntD(vxp,vyp); % 0
           if ((vxp+m1_xpp) <= RESX) & ((vyp+m1_ypp) <= RESY) & (vxp+m1_xpp)>= 1 & (vyp+m1_ypp)>= 1
           ntD2(uint16(vxp+m1_xpp),uint16(vyp+m1_ypp))=ntD(vxp,vyp); % 1
       end
           if ((vxp+m2_xpp) <= RESX) & ((vyp+m2_ypp) <= RESY) & (vxp+m2_xpp)>= 1 & (vyp+m2_ypp)>= 1
           ntD2(uint16(vxp+m2_xpp),uint16(vyp+m2_ypp))=ntD(vxp,vyp); % 2
       end
           if ((vxp+m3_xpp) <= RESX) & ((vyp+m3_ypp) <= RESY) & (vxp+m3_xpp)>= 1 & (vyp+m3_ypp)>= 1
           ntD2(uint16(vxp+m3_xpp),uint16(vyp+m3_ypp))=ntD(vxp,vyp); % 3
       end
           if ((vxp-m1_xpp) <= RESX) & ((vyp-m1_ypp) <= RESY) & (vxp-m1_xpp)>= 1 & (vyp-m1_ypp)>= 1
           ntD2(uint16(vxp-m1_xpp),uint16(vyp-m1_ypp))=ntD(vxp,vyp); % -1
       end
           if ((vxp-m2_xpp) <= RESX) & ((vyp-m2_ypp) <= RESY) & (vxp-m2_xpp)>= 1 & (vyp-m2_ypp)>= 1
           ntD2(uint16(vxp-m2_xpp),uint16(vyp-m2_ypp))=ntD(vxp,vyp); % -2
       end
           if (vxp-m3_xpp) <= RESX & (vyp-m3_ypp) <= RESY & (vxp-m3_xpp)>= 1 & (vyp-m3_ypp)>= 1
           ntD2(uint16(vxp-m3_xpp),uint16(vyp-m3_ypp))=ntD(vxp,vyp); % -3        
       end
       end
    end
end
imwrite(ntD2,'CGH3_ntD2.bmp','bmp');
% 180 degree rotation
for vxp = 1:RESX
    for vyp = 1:RESY
        if(ntD2(vxp,vyp) >= 1)
           ntD3(vxp,vyp)=ntD2(vxp,vyp); % 0
           if ((vxp+m1_xppp) <= RESX) & ((vyp+m1_yppp) <= RESY) & (vxp+m1_xppp)>= 1 & (vyp+m1_yppp)>= 1
           ntD3(uint16(vxp+m1_xppp),uint16(vyp+m1_yppp))=ntD2(vxp,vyp); % 1
       end
           if ((vxp+m2_xppp) <= RESX) & ((vyp+m2_yppp) <= RESY) & (vxp+m2_xppp)>= 1 & (vyp+m2_yppp)>= 1
           ntD3(uint16(vxp+m2_xppp),uint16(vyp+m2_yppp))=ntD2(vxp,vyp); % 2
       end
           if ((vxp+m3_xppp) <= RESX) & ((vyp+m3_yppp) <= RESY) & (vxp+m3_xppp)>= 1 & (vyp+m3_yppp)>= 1
           ntD3(uint16(vxp+m3_xppp),uint16(vyp+m3_yppp))=ntD2(vxp,vyp); % 3
       end
           if ((vxp-m1_xppp) <= RESX) & ((vyp-m1_yppp) <= RESY) & (vxp-m1_xppp)>= 1 & (vyp-m1_yppp)>= 1
           ntD3(uint16(vxp-m1_xppp),uint16(vyp-m1_yppp))=ntD2(vxp,vyp); % -1
       end
           if ((vxp-m2_xppp) <= RESX) & ((vyp-m2_yppp) <= RESY) & (vxp-m2_xppp)>= 1 & (vyp-m2_yppp)>= 1
           ntD3(uint16(vxp-m2_xppp),uint16(vyp-m2_yppp))=ntD2(vxp,vyp); % -2
       end
           if ((vxp-m3_xppp) <= RESX) & ((vyp-m3_yppp) <= RESY) & (vxp-m3_xppp)>= 1 & (vyp-m3_yppp)>= 1
           ntD3(uint16(vxp-m3_xppp),uint16(vyp-m3_yppp))=ntD2(vxp,vyp); % -3        
       end
       end
    end
end

end
end
imwrite(ntD3,'CGH3_ntD3.bmp','bmp');