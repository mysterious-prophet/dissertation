function XI=INTERPOL3D(X,dx,dy,dz,spacing)
% original size
[m,n,h]=size(X);
% original coordinates
x=(0:n-1)*dx;
y=(0:m-1)*dy;
z=(0:h-1)*dz;
% new dimensions diminished by one
N=floor(dx/spacing*(n-1));
M=floor(dy/spacing*(m-1));
H=floor(dz/spacing*(h-1));
% cubic interpolation grid of given positive spacing
[XQ,YQ,ZQ]=meshgrid(0:N,0:M,0:H);
% linear 3D interpolation
XI=interp3(x,y,z,X,XQ*spacing,YQ*spacing,ZQ*spacing,'linear');