function Cm = backproj(Fm);
% backproj(Fm)
% backproject the filtered image Fm
% Fm(p_i, theta_j) = g(omega_j, p_i) 
%   where omega_j = (cos(theta_j),sin(theta_j))
RES = length(Fm);
mins = -1.0;
maxs = 1.0;
step = (maxs - mins)/RES;
p = sqrt(2.0)*[mins+step/2:step:maxs-step/2];
theta = [2*pi/RES:2*pi/RES:2*pi];
for i=1:RES
    for j=1:RES
        x  = mins + i*step - step/2;
        y  = mins + j*step - step/2;        
        pw = x*cos(theta) + y*sin(theta);                             
        g = interp1(p,Fm',pw);
        Cm(i,j) = (2*pi/RES)*trace(g);
    end
end    
