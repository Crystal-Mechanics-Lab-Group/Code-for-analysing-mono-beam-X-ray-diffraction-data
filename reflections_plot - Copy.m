%storing variables for given combination at a given depth
B_depth = zeros(3,3,100,100);
U_depth = zeros(3,3,100,100);
T_depth = zeros(3,3,100,100);
A_depth = zeros(3,3,100,100);
strain_depth = zeros(3,3,100,100);
vol_strain_depth = zeros(3,3,100,100);
dev_strain_depth = zeros(3,3,100,100);

%surface of peaks at reconstructed depth as followss:
d471 = 12;
d435 = 7;
d354 = 9;
d804 = 9;
d635 = 10;
d554 = 10;

num_of_refl = 1:6;

surf_471 = -4.5;
surf_435 = -7;
surf_354 = -6;
surf_804 = -6;
surf_635 = -5;
surf_554 = -5;

%%
cd('C:\Users\sedm5079\Documents\Year2\W1Re3000He300Cbg\using 3d gauss fit');
qwidth_354 = dlmread('Q_Positions_Qwidth_354.txt');
qwidth_435 = dlmread('Q_Positions_Qwidth_435.txt');
qwidth_471 = dlmread('Q_Positions_Qwidth_471.txt');
qwidth_804 = dlmread('Q_Positions_Qwidth_804.txt');
qwidth_635 = dlmread('Q_Positions_Qwidth_635.txt');
qwidth_554 = dlmread('Q_Positions_Qwidth_554.txt');


for i=1:size(qwidth_354,1)
    if qwidth_354(i,1) == surf_354
        pos_354 = i;
        break;
    end
end

for i=1:size(qwidth_435,1)
    if qwidth_435(i,1) == surf_435
        pos_435 = i;
        break;
    end
end

for i=1:size(qwidth_471,1)
    if qwidth_471(i,1) == surf_471
        pos_471 = i;
        break;
    end
end

for i=1:size(qwidth_804,1)
    if qwidth_804(i,1) == surf_804
        pos_804 = i;
        break;
    end
end

for i=1:size(qwidth_635,1)
    if qwidth_635(i,1) == surf_635
        pos_635 = i;
        break;
    end
end

for i=1:size(qwidth_554,1)
    if qwidth_554(i,1) == surf_554
        pos_554 = i;
        break;
    end
end


%%

    
%reflections at depth
for depth=0

    %peak354
    str1 = 'C:\Users\sedm5079\Documents\Year2\W1Re3000He300Cbg\peak354\depth';
    i= d354 + depth;
    s = strcat(str1,num2str(i));
    cd(s);
    Q3D = dlmread('Qspace3D.txt');
    QZ_coord = dlmread('Qspace3DCorners.txt');
    cd('C:\Users\sedm5079\Documents\Year2\W1Re3000He300Cbg\using 3d gauss fit');
    width = 4*(qwidth_354(pos_354+depth,2));
    if width ==0
        width = 1e-3;
    end
    power = 0.5;
    qspace_Centre354_gauss(depth+1,:) = qspace_centre(Q3D,QZ_coord,width,power);


   
    %peak435
    str1 = 'C:\Users\sedm5079\Documents\Year2\W1Re3000He300Cbg\peak435\depth';
    i= d435 + depth;
    s = strcat(str1,num2str(i));
    cd (s);
    Q3D = dlmread('Qspace3D.txt');
    QZ_coord = dlmread('Qspace3DCorners.txt');
    cd('C:\Users\sedm5079\Documents\Year2\W1Re3000He300Cbg\using 3d gauss fit');
    width = 4*(qwidth_435(pos_435+depth,2));
    if width ==0
        width = 1e-3;
    end
    power = 0.5;
    qspace_Centre435_gauss(depth+1,:) = qspace_centre(Q3D,QZ_coord,width,power);


    %peak471
    str1 = 'C:\Users\sedm5079\Documents\Year2\W1Re3000He300Cbg\peak471\depth';
    i= d471 + depth;
    s = strcat(str1,num2str(i));
    cd (s);
    Q3D = dlmread('Qspace3D.txt');
    QZ_coord = dlmread('Qspace3DCorners.txt');
    cd('C:\Users\sedm5079\Documents\Year2\W1Re3000He300Cbg\using 3d gauss fit');
    width = 4*(qwidth_471(pos_471+depth,2));
    if width ==0
        width = 1e-3;
    end
    power = 0.5;
    qspace_Centre471_gauss(depth+1,:) = qspace_centre(Q3D,QZ_coord,width,power);


    %peak804
    str1 = 'C:\Users\sedm5079\Documents\Year2\W1Re3000He300Cbg\peak804\depth';
    i= d804 + depth;
    s = strcat(str1,num2str(i));
    cd (s);
    Q3D = dlmread('Qspace3D.txt');
    QZ_coord = dlmread('Qspace3DCorners.txt');
    cd('C:\Users\sedm5079\Documents\Year2\W1Re3000He300Cbg\using 3d gauss fit');
    width = 4*(qwidth_804(pos_804+depth,2));
    if width ==0
        width = 1e-3;
    end
    power = 0.5;   %4
    qspace_Centre804_gauss(depth+1,:) = qspace_centre(Q3D,QZ_coord,width,power);
  
    %peak635
    str1 = 'C:\Users\sedm5079\Documents\Year2\W1Re3000He300Cbg\peak635\depth';
    i= d635 + depth;
    s = strcat(str1,num2str(i));
    cd (s);
    Q3D = dlmread('Qspace3D.txt');
    QZ_coord = dlmread('Qspace3DCorners.txt');
    cd('C:\Users\sedm5079\Documents\Year2\W1Re3000He300Cbg\using 3d gauss fit');
    width = 4*(qwidth_635(pos_635+depth,2));
    if width ==0
        width = 1e-3;
    end
    power = 0.5;   %4
    qspace_Centre635_gauss(depth+1,:) = qspace_centre(Q3D,QZ_coord,width,power);
    
    %peak554
    str1 = 'C:\Users\sedm5079\Documents\Year2\W1Re3000He300Cbg\peak554\depth';
    i= d554 + depth;
    s = strcat(str1,num2str(i));
    cd (s);
    Q3D = dlmread('Qspace3D.txt');
    QZ_coord = dlmread('Qspace3DCorners.txt');
    cd('C:\Users\sedm5079\Documents\Year2\W1Re3000He300Cbg\using 3d gauss fit');
    width = 4*(qwidth_554(pos_804+depth,2));
    if width ==0
        width = 1e-3;
    end
    power = 0.5;   %4
    qspace_Centre554_gauss(depth+1,:) = qspace_centre(Q3D,QZ_coord,width,power);
    
    all_refl_qcenter (1:3,1) = qspace_Centre354_gauss(depth+1,1:3);
    all_refl_qcenter (1:3,2) = qspace_Centre435_gauss(depth+1,1:3);
    all_refl_qcenter (1:3,3) = qspace_Centre471_gauss(depth+1,1:3);
    all_refl_qcenter (1:3,4) = qspace_Centre804_gauss(depth+1,1:3);
    all_refl_qcenter (1:3,5) = qspace_Centre635_gauss(depth+1,1:3);
    all_refl_qcenter (1:3,6) = qspace_Centre554_gauss(depth+1,1:3);

    
    all_refl_hkl (1:3,1) = [3 5 4];
    all_refl_hkl (1:3,2) = [4 3 5];
    all_refl_hkl (1:3,3) = [4 7 1];
    all_refl_hkl (1:3,4) = [8 0 4];       
    all_refl_hkl (1:3,5) = [6 3 5];       
    all_refl_hkl (1:3,6) = [5 5 4];       

    
    cd('C:\Users\sedm5079\Documents\Year2\W1Re3000He300Cbg');
%     UB*all_refl_hkl = all_refl_qcenter OR UB =
    UB = all_refl_qcenter/all_refl_hkl;
    [U, B] = UB2UandB_2(UB);
    [A, A0] = AfromB(B);
    
    %transformation matrix
    T = A/A0;
    I = [1 0 0; 0 1 0; 0 0 1];
    strain = 0.5*(T+transpose(T))-I;
    i = (1/3)*trace(strain);
    vol_strain = [i 0 0; 0 i 0; 0 0 i];
    dev_strain = strain - vol_strain;
    
    for i=depth+1
        B_depth_wo_comb_gauss(:,:,i) = squeeze(B);
        U_depth_wo_comb_gauss(:,:,i) = squeeze(U);
        T_depth_wo_comb_gauss(:,:,i) = squeeze(T);
        A_depth_wo_comb_gauss(:,:,i) = squeeze(A);
        strain_depth_wo_comb_gauss(:,:,i) = squeeze(strain);
        vol_strain_depth_wo_comb_gauss(:,:,i) = squeeze(vol_strain);
        dev_strain_depth_wo_comb_gauss(:,:,i) = squeeze(dev_strain);
    end
    
   
end



%%
% rotating by U

for i=1:31
    rot = rotx(135)*U_depth_wo_comb_gauss(:,:,i);
    strain_depth_rot(:,:,i) = rot*strain_depth_wo_comb_gauss(:,:,i)*transpose(rot);
    vol_strain_depth_rot(:,:,i) = rot*vol_strain_depth_wo_comb_gauss(:,:,i)*transpose(rot);
    dev_strain_depth_rot(:,:,i) = rot*dev_strain_depth_wo_comb_gauss(:,:,i)*transpose(rot);
end

%%
%Considering the average of values below 7um depth to be 0, because it is
%the substrate. So every line is individually being normalized to this
%average value
for i=1:3
    strain_depth_rot(i,i,1:31) = strain_depth_rot(i,i,1:31)+(0-mean(strain_depth_rot(i,i,23:31)));
%     strain_depth_rot(i,i,14:27)=0;
end

strain_depth_rot(1,2,1:31) = strain_depth_rot(1,2,1:31)+(0-mean(strain_depth_rot(1,2,23:31)));
strain_depth_rot(1,3,1:31) = strain_depth_rot(1,3,1:31)+(0-mean(strain_depth_rot(1,3,23:31)));
strain_depth_rot(2,3,1:31) = strain_depth_rot(2,3,1:31)+(0-mean(strain_depth_rot(2,3,23:31)));


%%
%plotting principle strains 
figure
depth = linspace(0,15,31);
depth = depth*cosd(45);
lim = 2e-3;

    for i=1:3
        clear strain;
        strain(1,:) = squeeze(strain_depth_rot(i,i,1:31));
        p1 = plot(depth,strain);
        ylim([-lim,lim]);
        xlim([0,10]);
        p1.LineWidth = 2;
        hold on;
        xlabel('depth (\mum)','fontsize', 20,'FontWeight','bold');
        ylabel('strain','fontsize', 20,'FontWeight','bold');
    end
    h1 = legend('\epsilon_{xx}','\epsilon_{yy}','\epsilon_{zz}');
    set(h1,'fontsize',18);
    set(gca,'fontsize',18);
    title(sprintf('Principle Strain'),'fontsize',15);

%%
%plotting shear strains
% use strain_depth_rot
figure
clear strain;
depth = linspace(0,15,31);
depth = depth*cosd(45);
clear strain;
strain(1,:) = squeeze(strain_depth_rot(1,2,1:31));
ylim([-lim,lim]);
xlim([0,10]);
p1 = plot(depth,strain);
p1.LineWidth = 2;
hold on;
xlabel('depth (\mum)','fontsize', 20,'FontWeight','bold');
ylabel('strain','fontsize', 20,'FontWeight','bold');
clear strain;
strain(1,:) = squeeze(strain_depth_rot(1,3,1:31));
p2 = plot(depth,strain);
ylim([-lim,lim]);
xlim([0,10]);
p2.LineWidth = 2;
hold on;
clear strain;
strain(1,:) = squeeze(strain_depth_rot(2,3,1:31));
p3 = plot(depth,strain);
ylim([-lim,lim]);
xlim([0,10]);
p3.LineWidth = 2;
h1 = legend('\epsilon_{xy}','\epsilon_{xz}','\epsilon_{yz}');
set(h1,'fontsize',13);
set(gca,'fontsize',18);
title(sprintf('Shear Strain'),'fontsize',15);




