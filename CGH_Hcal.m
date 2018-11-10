function FPA = CGH_Hcal(G,part,FPA_size)
try
% 5 AUG 04 by L. Sethaphong
% forward projection along different axes defined by theta
% 3D to 2D projection from phi, theta
% let phi be 0 and pi/4
% let the distance to the FPA be 100 volume elements
% let the lamda volume be 10 units in height
% ===================================================
RES = length(G);
npi = 3.14159265358932384626433;
max = FPA_size;
phi = npi/6; 
nF = zeros(max,max,3); % create a new tiff file
%imwrite(ntF,fname,'tiff');
%nF = imread(fname); % read back the empty tiff file
vHeight = RES*2;
lumens = 0; % intensity scaling
% main projections
middle = max/2 + 1;
origin_D = RES/2 + 1;
origin_X = middle;
origin_Y = middle;
partitions = 2*part;
distMap = zeros(RES,RES);
distMap = rgb2wave(G);

% center mapping
    for m=1:RES
        new_X = uint16(origin_X + (m-origin_D));
        for n=1:RES
            new_Y = uint16(origin_Y + (n-origin_D));
            for k=1:3
                nF(new_X,new_Y,k) = uint8(0.85*double(G(m,n,k)));
            end
            %            nF(new_X,new_Y,k) = G(m,n,k);
        end
    end
    
% CGH SIM Projections
    for s=1:2
        lamda = vHeight*tan(phi*s);
%        separation = 15*(4-k)*sin(phi*s) + lamda;
        for theta=1:(partitions*s)
             rot_theta = ((theta-1)*npi/(part*s));
%             x_sep = separation*cos(rot_theta);
%             y_sep = separation*sin(rot_theta);
             x_sep = cos(rot_theta);
             y_sep = sin(rot_theta);
            for i=1:RES
                for j=1:RES
                    % was 2*distMap
                    separation = lamda + sin(phi*s)*distMap(i,j);
                    new_X = uint16(origin_X + (i-origin_D) + x_sep*separation);
                    new_Y = uint16(origin_Y + (j-origin_D) + y_sep*separation);
                    % need scaling
                    if new_X < 1
                        new_X = 1;
                    elseif new_Y < 1
                        new_Y = 1;
                    end
                    if new_X < length(nF) & new_Y < length(nF)
                       for k=1:3
                           if(double(G(i,j,k))+ double(nF(new_X,new_Y,k))) > 254
                            nF(new_X,new_Y,k) = uint8(255);
                           else
                            nF(new_X,new_Y,k) = uint8(0.5*cos(phi)*double(G(i,j,k))+0.5*cos(phi)*double(nF(new_X,new_Y,k)));
                           end
                        end
                    end
                end
            end
        end
    end
% final output
FPA=nF;
%imwrite(nF,fname,'tiff');
catch
    i
    j
    k
    new_X
    new_Y
end