function[R,Xp,Yp]=radon3(I,theta,phi)
% Radon3 Compute 3-D Radon Transform
% ----------------------------------
% v 1.0b 3/09/1997
% by P. Cipollini & P. Challenor
[yi,xi,zi] = size(I);
pi = acos(-1);
x=-floor(xi/2):floor(xi/2);
y=-floor(yi/2):floor(yi/2);
z=-floor(zi/2):floor(zi/2);
[X Y Z]=meshgrid(x,y,z);
theta=(theta/180)*pi; 
phi=(phi/180)*pi;
X1=cos(phi)*X+sin(phi)*Y;
Y1=-cos(theta)*sin(phi)*X+cos(theta)*cos(phi)*Y+sin(theta)*Z;
Xr=round(X1(:));
Yr=round(Y1(:));
my=min(Yr)-1;
mx=min(Xr)-1;
R=zeros(max(Yr)-my,max(Xr)-mx);
for i=1:length(I(:)),
	R(Yr(i)-my,Xr(i)-mx)=R(Yr(i)-my,Xr(i)-mx)+I(i);
end
Xp=min(Xr):max(Xr);
Yp=min(Yr):max(Yr);