function y=gauss3D(X,x)
% This function is called by leasqr
% x is a vector which contains the coefficients of the
% equation. X and P are the option data sets that were
% passed to leasqr. The first column in X contains xcoords, the second
% ycoords

%edited by F Hofmann 27th August 2014

%inputs:
% A-peak height  B-centre (B1, B2 and B3)  C-width (C)  D-background


B1 = x(1); %center in dimension 1
B2 = x(2); %center in dimension 2
B3 = x(3); %center in dimension 3
C = x(4); %width
A = x(5); %prefactor
D = x(6); %background

%coordinates:
coor1 = X(:,1);
coor2 = X(:,2);
coor3 = X(:,3);

%compute 3D Gaussian:

y = A.*...
    exp(-((coor1-B1).^2)./(2*C^2)).*...
    exp(-((coor2-B2).^2)./(2*C^2)).*...
    exp(-((coor3-B3).^2)./(2*C^2)).*...
    + D;

