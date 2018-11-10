function Im = radon(N, C, R, D, RES)
% radon(N, C, R, D, RES)
% analytic Radon transform of N overlapping circles
% N = number of circles
% C = centers of the circles (x_i,y_i) = (C(i,1),C(i,2))
% R = radii of the circles R(i)
% D = density of the circles D(i)
% RES = resolution
% ------------------------------------------------------------
mins = -1.0;
maxs =  1.0;
step = (maxs - mins)/RES;
% theta and p
Im = zeros(RES,RES);
theta = [2*pi/RES:2*pi/RES:2*pi];
pp  = sqrt(2.0)*[step/2:step:maxs-step/2];
pn =  sqrt(2.0)*[mins+step/2:step:-step/2];
% for each row k
for k=1:RES/2
    % sum over all the overlapping circles
    for i=1:N
        tmp0 =  abs( pn(k)*ones(1,RES) - C(i,1)*cos(theta) - C(i,2)*sin(theta) );
        tmp1 =  R(i)^2*ones(1,RES) - tmp0.^2;
        % get rid of the negative terms
        tmp2 =  (tmp1 > zeros(1,RES)).* tmp1;
        Im(:,k) = Im(:,k) + (2*D(i)*sqrt(tmp2))';
    end    
end   
for k=1:RES/2
    % sum over all the overlapping circles
    for i=1:N
        tmp0 =  abs( pp(k)*ones(1,RES) - C(i,1)*cos(theta) - C(i,2)*sin(theta) );
        tmp1 =  R(i)^2*ones(1,RES) - tmp0.^2;
        % get rid of the negative terms
        tmp2 =  (tmp1 > zeros(1,RES)).* tmp1;
        Im(:,k+RES/2) = Im(:,k+RES/2) + (2*D(i)*sqrt(tmp2))';
    end    
end  
