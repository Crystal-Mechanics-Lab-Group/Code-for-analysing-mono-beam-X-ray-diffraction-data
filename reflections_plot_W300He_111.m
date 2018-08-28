%storing variables for given combination at a given depth
B_depth = zeros(3,3,100,100);
U_depth = zeros(3,3,100,100);
T_depth = zeros(3,3,100,100);
A_depth = zeros(3,3,100,100);
strain_depth = zeros(3,3,100,100);
vol_strain_depth = zeros(3,3,100,100);
dev_strain_depth = zeros(3,3,100,100);

%surface of peaks at reconstructed depth as followss:
d075 = 20;
d770 = 15;
d78m1 = 16;
d345 = 17;
d444 = 16;
d653 = 16;
d464 = 17;

surf_075 = -0.5;
surf_345 = -2;
surf_444 = -2.5;
surf_464 = -2;
surf_653 = -2.5;
surf_770 = -3;
surf_78m1 = -2.5;


%%
cd('W300He_111');
qwidth_075 = dlmread('Qwidth075.txt');
qwidth_770 = dlmread('Qwidth770.txt');
qwidth_78m1 = dlmread('Qwidth78m1.txt');
qwidth_345 = dlmread('Qwidth345.txt');
qwidth_444 = dlmread('Qwidth444.txt');
qwidth_653 = dlmread('Qwidth653.txt');
qwidth_464 = dlmread('Qwidth464.txt');


%%
for i=1:size(qwidth_075,1)
    if qwidth_075(i,1) == surf_075
        pos_075 = i;
        break;
    end
end

for i=1:size(qwidth_770,1)
    if qwidth_770(i,1) == surf_770
        pos_770 = i;
        break;
    end
end

for i=1:size(qwidth_78m1,1)
    if qwidth_78m1(i,1) == surf_78m1
        pos_78m1 = i;
        break;
    end
end

for i=1:size(qwidth_345,1)
    if qwidth_345(i,1) == surf_345
        pos_345 = i;
        break;
    end
end

for i=1:size(qwidth_444,1)
    if qwidth_444(i,1) == surf_444
        pos_444 = i;
        break;
    end
end

for i=1:size(qwidth_653,1)
    if qwidth_653(i,1) == surf_653
        pos_653 = i;
        break;
    end
end

for i=1:size(qwidth_464,1)
    if qwidth_464(i,1) == surf_464
        pos_464 = i;
        break;
    end
end

cd('..\');

%%
    
%reflections at depth
for depth=0:30

    %peak075
    str1 = 'W300He_111\peak075_pos2\depth';
    i= d075 + depth;
    s = strcat(str1,num2str(i));
    cd(s);
    Q3D = dlmread('Qspace3D.txt');
    QZ_coord = dlmread('Qspace3DCorners.txt');
    cd('..\..\..\');
    width = 4*(qwidth_075(pos_075+depth,2));
    if width ==0
        width = 1e-3;
    end
    power = 0.5;
    qspace_Centre075_gauss(depth+1,:) = qspace_centre(Q3D,QZ_coord,width,power);


   
    %peak770
    str1 = 'W300He_111\peak770_pos2\depth';
    i= d770 + depth;
    s = strcat(str1,num2str(i));
    cd (s);
    Q3D = dlmread('Qspace3D.txt');
    QZ_coord = dlmread('Qspace3DCorners.txt');
    cd('..\..\..\');
    width = 4*(qwidth_770(pos_770+depth,2));
    if width ==0
        width = 1e-3;
    end
    power = 0.5;
    qspace_Centre770_gauss(depth+1,:) = qspace_centre(Q3D,QZ_coord,width,power);


    %peak78-1
    str1 = 'W300He_111\peak78-1\depth';
    i= d78m1 + depth;
    s = strcat(str1,num2str(i));
    cd (s);
    Q3D = dlmread('Qspace3D.txt');
    QZ_coord = dlmread('Qspace3DCorners.txt');
    cd('..\..\..\');
    width = 4*(qwidth_78m1(pos_78m1+depth,2));
    if width ==0
        width = 1e-3;
    end
    power = 0.5;
    qspace_Centre78m1_gauss(depth+1,:) = qspace_centre(Q3D,QZ_coord,width,power);


    %peak345
    str1 = 'W300He_111\peak345\depth';
    i= d345 + depth;
    s = strcat(str1,num2str(i));
    cd (s);
    Q3D = dlmread('Qspace3D.txt');
    QZ_coord = dlmread('Qspace3DCorners.txt');
    cd('..\..\..\');
    width = 4*(qwidth_345(pos_345+depth,2));
    if width ==0
        width = 1e-3;
    end
    power = 0.5;   %4
    qspace_Centre345_gauss(depth+1,:) = qspace_centre(Q3D,QZ_coord,width,power);


    %peak444
    str1 = 'W300He_111\peak444_pos2\depth';
    i= d444 + depth;
    s = strcat(str1,num2str(i));
    cd (s);
    Q3D = dlmread('Qspace3D.txt');
    QZ_coord = dlmread('Qspace3DCorners.txt');
    cd('..\..\..\');
    width = 4*(qwidth_444(pos_444+depth,2));
    if width ==0
        width = 1e-3;
    end
    power = 0.5;%8;%5
    qspace_Centre444_gauss(depth+1,:) = qspace_centre(Q3D,QZ_coord,width,power);


    %peak653
    str1 = 'W300He_111\peak653_pos2\depth';
    i= d653 + depth;
    s = strcat(str1,num2str(i));
    cd (s);
    Q3D = dlmread('Qspace3D.txt');
    QZ_coord = dlmread('Qspace3DCorners.txt');
    cd('..\..\..\');
    width = 4*(qwidth_653(pos_653+depth,2));
    if width ==0
        width = 1e-3;
    end
    power = 0.5;%12;%5
    qspace_Centre653_gauss(depth+1,:) = qspace_centre(Q3D,QZ_coord,width,power);
     
    
    
    %peak464
    str1 = 'W300He_111\peak464_pos2\depth';
    i= d464 + depth;
    s = strcat(str1,num2str(i));
    cd (s);
    Q3D = dlmread('Qspace3D.txt');
    QZ_coord = dlmread('Qspace3DCorners.txt');
    cd('..\..\..\');
    width = 4*(qwidth_464(pos_464+depth,2));
    if width ==0
        width = 1e-3;
    end
    power = 0.5;%8;%5
    qspace_Centre464_gauss(depth+1,:) = qspace_centre(Q3D,QZ_coord,width,power);
  
    
    all_refl_qcenter (1:3,1) = qspace_Centre075_gauss(depth+1,1:3);
    all_refl_qcenter (1:3,2) = qspace_Centre345_gauss(depth+1,1:3);
    all_refl_qcenter (1:3,3) = qspace_Centre770_gauss(depth+1,1:3);
    all_refl_qcenter (1:3,4) = qspace_Centre78m1_gauss(depth+1,1:3);
    all_refl_qcenter (1:3,5) = qspace_Centre653_gauss(depth+1,1:3);
    all_refl_qcenter (1:3,6) = qspace_Centre464_gauss(depth+1,1:3);
    all_refl_qcenter (1:3,7) = qspace_Centre444_gauss(depth+1,1:3);
    
    all_refl_hkl (1:3,1) = [0 7 5];
    all_refl_hkl (1:3,2) = [3 4 5];
    all_refl_hkl (1:3,3) = [7 7 0];
    all_refl_hkl (1:3,4) = [7 8 -1];       
    all_refl_hkl (1:3,5) = [6 5 3];
    all_refl_hkl (1:3,6) = [4 6 4];
    all_refl_hkl (1:3,7) = [4 4 4];
    
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
%rotating by U

for i=1:31
    rot = rotx(135)*U_depth_wo_comb_gauss(:,:,i);
    strain_depth_rot(:,:,i) = rot*strain_depth_wo_comb_gauss(:,:,i)*transpose(rot);
    vol_strain_depth_rot(:,:,i) = rot*vol_strain_depth_wo_comb_gauss(:,:,i)*transpose(rot);
    dev_strain_depth_rot(:,:,i) = rot*dev_strain_depth_wo_comb_gauss(:,:,i)*transpose(rot);
end

%%
%Considering the average of values below 7 um depth to be 0, because it is
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
lim = 1e-3;
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
set(h1,'fontsize',13);
set(gca,'fontsize',18);
title(sprintf('Principle Strain'),'fontsize',15);

%%
%plotting shear strains
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



