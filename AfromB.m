function [A, A0] = AfromB(B)
%%
% %magnitude of unit vectors in reciprocal lattice
B1 = B(:,1);
B2 = B(:,2);
B3 = B(:,3);

%volume of unit cell
V = dot(B1,cross(B2,B3));

%unit vectors in real space
A(:,1) = transpose((2*pi)*(cross(B2,B3)./V));
A(:,2) = transpose((2*pi)*(cross(B3,B1)./V));
A(:,3) = transpose((2*pi)*(cross(B1,B2)./V));
A = transpose(A);

%contructing A matrix for unstrained cell
a = 0.31652;
AA1 = [a 0 0];
AA2 = [0 a 0];
AA3 = [0 0 a];
V = dot(AA1, cross(AA2,AA3));
b = cross(AA1,AA2)/V;
angle = pi/2;
A0(1,1) = a;
% A0(1,2) = a*cos(angle);
A0(1,2) = 0;
% A0(1,3) = a*cos(angle);
A0(1,3) = 0;
A0(2,1) = 0;
A0(2,2) = a*sin(angle);
% A0(2,3) = -a*sin(angle)*cos(angle);
A0(2,3) = 0;
A0(3,1) = 0;
A0(3,2) = 0;
A0(3,3) = 1/b(1,3);

%%
%IGNORE
% %magnitude of unit vectors in reciprocal lattice
% b1 = B(1,1);
% %beta3 = arccos(B1.B2/abs(B1)*abs(B2)) so B1 = [b1 0 0] and B2 = [0 b2 0]
% %so dot(B1,B2) will be 0. So beta3 = acos(0). Same applies for beta2 and
% %beta1
% beta3 = abs(atan(B(2,2)/B(1,2)));
% b2 = B(2,2)/sin(beta3);
% b3 = B(1,3)/sin(beta3); %check
% 
% %unit vectors in reciprocal lattice
% % B1 = [b1 0 0];
% % B2 = [0 b2 0];
% % B3 = [0 0 b3];

%magnitude of unit vectors in real space
% a1 = A1(1,1);
% a2 = A2(1,2);
% a3 = A3(1,3);
% 
% %constructing A matrix in real space
% %a1 = (2*pi)/(cross(b2,b3)/V). Now V = dot(b1,cross(b2,b3))
% 
% A(1,1) = a1;
% % A(1,2) = a2*cos(beta3);
% A(1,2) = 0;
% % A(1,3) = a3*cos(beta3);
% A(1,3) = 0;
% A(2,1) = 0;
% A(2,2) = a2*sin(beta3);
% % A(2,3) = -a3*sin(beta3)*cos(beta3);
% A(2,3) = 0;
% A(3,1) = 0;
% A(3,2) = 0;
% A(3,3) = a2; %because b3 comes out to be very strange when dividing B(1,3) by cos(beta3) because its cos 90 which is near 0

end

