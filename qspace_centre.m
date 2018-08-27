function qspace_centre = qspace_centre(Q3D,QZ_coord,width,power)

clear q_matrix;
clear q_matrix_new;
clear q3d_data;
clear qx;
clear qy;
clear qz;
clear qspace_centre;

q3d_data = Q3D(2:end,2:end);

%find qx
qx_start=Q3D(2,1);
for i=3:size(Q3D,1);
    if Q3D(i,1) <= qx_start
        break;
    end
end
qx = Q3D(2:(i-1),1);

%find qy
qy = transpose(Q3D(1,1:(size(Q3D,2)-1)));

%find qz
count=0;
for i=3:size(Q3D,1)
    if Q3D(i,1) == Q3D(2,1)
        count = count+1;
    end
end

stepsize = (QZ_coord(3,4) - QZ_coord(2,4))/(count-1);
qz_grid = meshgrid(QZ_coord(2,4):stepsize:QZ_coord(3,4));
qz = transpose(qz_grid(1,:));

%constructing matrix with qx in coulm 1, qy in col2, qz in col 3 and
%intensity in col 4
q3d_data(q3d_data==0) = NaN;
q3d_data(isnan(q3d_data)) = -5;

row=1;
q_matrix_row=1;
for i=1:size(qz,1)
    for j=1:size(qx,1)
        for k=1:size(qy,1)
            q_matrix(q_matrix_row,1) = qx(j,1); %storing x component
            q_matrix(q_matrix_row,2) = qy(k,1); %storing y component 
            q_matrix(q_matrix_row,3) = qz(i,1); %storing z component
            q_matrix(q_matrix_row,4) = q3d_data(row,k); %storing q_data component
            q_matrix_row = q_matrix_row +1;  
        end
        row=row+1;
    end
end

%removing all NaN or zero data from intensity
j=1;
for i=1:size(q_matrix,1)
    if q_matrix(i,4) ~= -5
        q_matrix_new(j,1) = q_matrix(i,1);
        q_matrix_new(j,2) = q_matrix(i,2);
        q_matrix_new(j,3) = q_matrix(i,3);
        q_matrix_new(j,4) = q_matrix(i,4);
        j=j+1;
    end
end

pos = q_matrix_new(:,1:3);
Idat = q_matrix_new(:,4);
 
% B1 = x(1); %center in dimension 1
% B2 = x(2); %center in dimension 2
% B3 = x(3); %center in dimension 3
% C = x(4); %width
% A = x(5); %prefactor
% D = x(6); %background
 
% make initial guess...
[Alt,k0]= max(Idat);
x0(1) = q_matrix_new(k0,1);
x0(2) = q_matrix_new(k0,2);
x0(3) = q_matrix_new(k0,3);
% x0(4) = guess for width -> rough exp width /10
x0(4) = width;
x0(5) = Alt;
x0(6) = 0;
w = Idat.^power;
cd('C:\Users\sedm5079\Documents\Year2\W300He300C_bg2\using 3d gauss fit');

[yfit,x,~,~,~,~]=leasqr(pos,Idat,x0,'gauss3D',0.0001,40,w);
qspace_centre = x(1:3,1)';

% figure
% scatter3(q_matrix_new(:,1),q_matrix_new(:,2),q_matrix_new(:,3),5,q_matrix_new(:,4));
% colorbar;
% xlabel('Qx')
% ylabel('Qy')
% zlabel('Qz')
% % % 
% qspace_centre = findcentre(q_matrix_new);
% hold on;
% scatter3(qspace_centre(1,1),qspace_centre(1,2), qspace_centre(1,3),20,'red','filled');
% scatter3(qspace_centre(1,1),qspace_centre(2,1), qspace_centre(3,1),20,'red','filled');

end


