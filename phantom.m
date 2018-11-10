function Pm = phantom(N, C, R, D, RES)
% phantom(N, C, R, D, RES)
% N = number of circles
% C = centers of the circles (x_i,y_i) = (C(i,1),C(i,2))
% R = radii of the circles R(i)
% D = density of the circles D(i)
% RES = resolution of the image
% generate a phantom images of N overlapping circle
% ==============================================================
Pm = zeros(RES,RES);
% for scaling purpose
mins = -1.0;
maxs = 1.0;
scale = RES/(maxs - mins);
% filling the desity of circles
for i=1:N
    xc = round((C(i,1)-mins)*scale);
    yc = round((C(i,2)-mins)*scale);
    rc = round(R(i)*scale);
    for x = (xc-rc):(xc+rc)
        for y = (yc-rc):(yc+rc)
            if ((x-xc)^2 + (y-yc)^2 < rc^2)
               Pm(x,y) = Pm(x,y) + D(i);    
            end   
        end
    end    
end      
    
    
