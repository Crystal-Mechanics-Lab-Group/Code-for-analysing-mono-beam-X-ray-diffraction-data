%storing variables for given combination at a given depth
B_depth = zeros(3,3,100,100);
U_depth = zeros(3,3,100,100);
T_depth = zeros(3,3,100,100);
A_depth = zeros(3,3,100,100);
strain_depth = zeros(3,3,100,100);
vol_strain_depth = zeros(3,3,100,100);
dev_strain_depth = zeros(3,3,100,100);

%surface of peaks at reconstructed depth as followss:
d5m29 = 16;
dm4m19 = 11;
d037 = 11;
d318 = 15;
d019 = 12;
d127 = 12;
d226 = 13;


surf_5m29 = -2;
surf_m4m19 = -4.5;
surf_037 = -4.5;
surf_318 = -2.5;
surf_019 = -4;
surf_127 = -4;
surf_226 = -3.5;

%%
cd('W1Re300He_001');
qwidth_5m29 = dlmread('Q_Positions_Qwidth_5-29.txt');
qwidth_m4m19 = dlmread('Q_Positions_Qwidth_-4-19.txt');
qwidth_037 = dlmread('Q_Positions_Qwidth_037.txt');
qwidth_318 = dlmread('Q_Positions_Qwidth_318.txt');
qwidth_019 = dlmread('Q_Positions_Qwidth_019.txt');
qwidth_127 = dlmread('Q_Positions_Qwidth_127.txt');
qwidth_226 = dlmread('Q_Positions_Qwidth_226.txt');


for i=1:size(qwidth_5m29,1)
    if qwidth_5m29(i,1) == surf_5m29
        pos_5m29 = i;
        break;
    end
end

for i=1:size(qwidth_m4m19,1)
    if qwidth_m4m19(i,1) == surf_m4m19
        pos_m4m19 = i;
        break;
    end
end

for i=1:size(qwidth_037,1)
    if qwidth_037(i,1) == surf_037
        pos_037 = i;
        break;
    end
end

for i=1:size(qwidth_318,1)
    if qwidth_318(i,1) == surf_318
        pos_318 = i;
        break;
    end
end

for i=1:size(qwidth_019,1)
    if qwidth_019(i,1) == surf_019
        pos_019 = i;
        break;
    end
end

for i=1:size(qwidth_127,1)
    if qwidth_127(i,1) == surf_127
        pos_127 = i;
        break;
    end
end

for i=1:size(qwidth_226,1)
    if qwidth_226(i,1) == surf_226
        pos_226 = i;
        break;
    end
end

cd('..\');

%%    
%reflections at depth
for depth=0:30

    %peak5-29
    str1 = 'W1Re300He_001\peak5-29\depth';
    i= d5m29 + depth;
    s = strcat(str1,num2str(i));
    cd(s);
    Q3D = dlmread('Qspace3D.txt');
    QZ_coord = dlmread('Qspace3DCorners.txt');
    cd('..\..\..\');
    width = 4*(qwidth_5m29(pos_5m29+depth,2));
    if width ==0
        width = 1e-3;
    end
    power = 0.5;
    qspace_Centre5m29_gauss(depth+1,:) = qspace_centre(Q3D,QZ_coord,width,power);

   
    %peak-4-19
    str1 = 'W1Re300He_001\peak-4-19\depth';
    i= dm4m19 + depth;
    s = strcat(str1,num2str(i));
    cd (s);
    Q3D = dlmread('Qspace3D.txt');
    QZ_coord = dlmread('Qspace3DCorners.txt');
    cd('..\..\..\');
    width = 4*(qwidth_m4m19(pos_m4m19+depth,2));
    if width ==0
        width = 1e-3;
    end
    power = 0.5;
    qspace_Centrem4m19_gauss(depth+1,:) = qspace_centre(Q3D,QZ_coord,width,power);


    %peak037
    str1 = 'W1Re300He_001\peak037\depth';
    i= d037 + depth;
    s = strcat(str1,num2str(i));
    cd (s);
    Q3D = dlmread('Qspace3D.txt');
    QZ_coord = dlmread('Qspace3DCorners.txt');
    cd('..\..\..\');
    width = 4*(qwidth_037(pos_037+depth,2));
    if width ==0
        width = 1e-3;
    end
    power = 0.5;
    qspace_Centre037_gauss(depth+1,:) = qspace_centre(Q3D,QZ_coord,width,power);


    %peak318
    str1 = 'W1Re300He_001\peak318\depth';
    i= d318 + depth;
    s = strcat(str1,num2str(i));
    cd (s);
    Q3D = dlmread('Qspace3D.txt');
    QZ_coord = dlmread('Qspace3DCorners.txt');
    cd('..\..\..\');
    width = 4*(qwidth_318(pos_318+depth,2));
    if width ==0
        width = 1e-3;
    end
    power = 0.5;   %4
    qspace_Centre318_gauss(depth+1,:) = qspace_centre(Q3D,QZ_coord,width,power);


    %peak019
    str1 = 'W1Re300He_001\peak019\depth';
    i= d019 + depth;
    s = strcat(str1,num2str(i));
    cd (s);
    Q3D = dlmread('Qspace3D.txt');
    QZ_coord = dlmread('Qspace3DCorners.txt');
    cd('..\..\..\');
    width = 4*(qwidth_019(pos_019+depth,2));
    if width ==0
        width = 1e-3;
    end
    power = 0.5;%8;%5
    qspace_Centre019_gauss(depth+1,:) = qspace_centre(Q3D,QZ_coord,width,power);
  
    %peak127
    str1 = 'W1Re300He_001\peak127\depth';
    i= d127 + depth;
    s = strcat(str1,num2str(i));
    cd (s);
    Q3D = dlmread('Qspace3D.txt');
    QZ_coord = dlmread('Qspace3DCorners.txt');
    cd('..\..\..\');
    width = 4*(qwidth_127(pos_127+depth,2));
    if width ==0
        width = 1e-3;
    end
    power = 0.5;%8;%5
    qspace_Centre127_gauss(depth+1,:) = qspace_centre(Q3D,QZ_coord,width,power);
    
    %peak226
    str1 = 'W1Re300He_001\peak226\depth';
    i= d226 + depth;
    s = strcat(str1,num2str(i));
    cd (s);
    Q3D = dlmread('Qspace3D.txt');
    QZ_coord = dlmread('Qspace3DCorners.txt');
    cd('..\..\..\');
    width = 4*(qwidth_226(pos_226+depth,2));
    if width ==0
        width = 1e-3;
    end
    power = 0.5;%8;%5
    qspace_Centre226_gauss(depth+1,:) = qspace_centre(Q3D,QZ_coord,width,power);
    
    all_refl_qcenter (1:3,1) = qspace_Centre5m29_gauss(depth+1,1:3);
    all_refl_qcenter (1:3,2) = qspace_Centrem4m19_gauss(depth+1,1:3);
    all_refl_qcenter (1:3,3) = qspace_Centre037_gauss(depth+1,1:3);
    all_refl_qcenter (1:3,4) = qspace_Centre318_gauss(depth+1,1:3);
    all_refl_qcenter (1:3,5) = qspace_Centre019_gauss(depth+1,1:3);
    all_refl_qcenter (1:3,6) = qspace_Centre127_gauss(depth+1,1:3);
    all_refl_qcenter (1:3,7) = qspace_Centre226_gauss(depth+1,1:3);

    
    all_refl_hkl (1:3,1) = [5 -2 9];
    all_refl_hkl (1:3,2) = [-4 -1 9];
    all_refl_hkl (1:3,3) = [0 3 7];
    all_refl_hkl (1:3,4) = [3 1 8];       
    all_refl_hkl (1:3,5) = [0 1 9];
    all_refl_hkl (1:3,6) = [1 2 7];
    all_refl_hkl (1:3,7) = [2 2 6];

    
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
%    
end

%%
%rotating by U
clear strain_depth_rot;

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
    strain_depth_rot(i,i,1:31) = strain_depth_rot(i,i,1:31)+(0-mean(strain_depth_rot(i,i,18:31)));
end

strain_depth_rot(1,2,1:31) = strain_depth_rot(1,2,1:31)+(0-mean(strain_depth_rot(1,2,18:31)));
strain_depth_rot(1,3,1:31) = strain_depth_rot(1,3,1:31)+(0-mean(strain_depth_rot(1,3,18:31)));
strain_depth_rot(2,3,1:31) = strain_depth_rot(2,3,1:31)+(0-mean(strain_depth_rot(2,3,18:31)));


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
depth = depth *cosd(45);
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



