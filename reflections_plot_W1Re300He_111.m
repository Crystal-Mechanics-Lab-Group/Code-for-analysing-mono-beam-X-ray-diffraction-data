%storing variables for given combination at a given depth
B_depth = zeros(3,3,100,100);
U_depth = zeros(3,3,100,100);
T_depth = zeros(3,3,100,100);
A_depth = zeros(3,3,100,100);
strain_depth = zeros(3,3,100,100);
vol_strain_depth = zeros(3,3,100,100);
dev_strain_depth = zeros(3,3,100,100);

%surface of peaks at reconstructed depth as followss:
d354 = 17;
d444 = 15;
d563 = 18;
d804 = 14;
d680 = 21;
d745 = 15;
d1024 = 17;

surf_354 = -1.5;
surf_444 = -2.5;
surf_563 = -1;
surf_804 = -3;
surf_680 = 0.5;
surf_745 = -2.5;
surf_1024 = -1.5;

%%
cd('W1Re300He_111');
qwidth_354 = dlmread('Qwidth354.txt');
qwidth_444 = dlmread('Qwidth444.txt');
qwidth_563 = dlmread('Qwidth563.txt');
qwidth_804 = dlmread('Qwidth804.txt');
qwidth_680 = dlmread('Qwidth680.txt');
qwidth_745 = dlmread('Qwidth745.txt');
qwidth_1024 = dlmread('Qwidth1024.txt');


for i=1:size(qwidth_354,1)
    if qwidth_354(i,1) == surf_354
        pos_354 = i;
        break;
    end
end

for i=1:size(qwidth_444,1)
    if qwidth_444(i,1) == surf_444
        pos_444 = i;
        break;
    end
end

for i=1:size(qwidth_563,1)
    if qwidth_563(i,1) == surf_563
        pos_563 = i;
        break;
    end
end

for i=1:size(qwidth_804,1)
    if qwidth_804(i,1) == surf_804
        pos_804 = i;
        break;
    end
end

for i=1:size(qwidth_680,1)
    if qwidth_680(i,1) == surf_680
        pos_680 = i;
        break;
    end
end

for i=1:size(qwidth_745,1)
    if qwidth_745(i,1) == surf_745
        pos_745 = i;
        break;
    end
end

for i=1:size(qwidth_1024,1)
    if qwidth_1024(i,1) == surf_1024
        pos_1024 = i;
        break;
    end
end

cd('..\');

%%
    
%reflections at depth
for depth=0:28

    %peak354
    str1 = 'W1Re300He_111\peak354\depth';
    i= d354 + depth;
    s = strcat(str1,num2str(i));
    cd(s);
    Q3D = dlmread('Qspace3D.txt');
    QZ_coord = dlmread('Qspace3DCorners.txt');
    cd('..\..\..\');
    width = 4*(qwidth_354(pos_354+depth,2));
    if width ==0
        width = 1e-3;
    end
    power = 0.5;
    qspace_Centre354_gauss(depth+1,:) = qspace_centre(Q3D,QZ_coord,width,power);

   
    %peak444
    str1 = 'W1Re300He_111\peak444\depth';
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
    power = 0.5;
    qspace_Centre444_gauss(depth+1,:) = qspace_centre(Q3D,QZ_coord,width,power);


    %peak563
    str1 = 'W1Re300He_111\peak563\depth';
    i= d563 + depth;
    s = strcat(str1,num2str(i));
    cd (s);
    Q3D = dlmread('Qspace3D.txt');
    QZ_coord = dlmread('Qspace3DCorners.txt');
    cd('..\..\..\');
    width = 4*(qwidth_563(pos_563+depth,2));
    if width ==0
        width = 1e-3;
    end
    power = 0.5;
    qspace_Centre563_gauss(depth+1,:) = qspace_centre(Q3D,QZ_coord,width,power);


    %peak804
    str1 = 'W1Re300He_111\peak804\depth';
    i= d804 + depth;
    s = strcat(str1,num2str(i));
    cd (s);
    Q3D = dlmread('Qspace3D.txt');
    QZ_coord = dlmread('Qspace3DCorners.txt');
    cd('..\..\..\');
    width = 4*(qwidth_804(pos_804+depth,2));
    if width ==0
        width = 1e-3;
    end
    power = 0.5;   %4
    qspace_Centre804_gauss(depth+1,:) = qspace_centre(Q3D,QZ_coord,width,power);


    %peak680
    str1 = 'W1Re300He_111\peak680\depth';
    i= d680 + depth;
    s = strcat(str1,num2str(i));
    cd (s);
    Q3D = dlmread('Qspace3D.txt');
    QZ_coord = dlmread('Qspace3DCorners.txt');
    cd('..\..\..\');
    width = 4*(qwidth_680(pos_680+depth,2));
    if width ==0
        width = 1e-3;
    end
    power = 0.5;%8;%5
    qspace_Centre680_gauss(depth+1,:) = qspace_centre(Q3D,QZ_coord,width,power);
  
    %peak745
    str1 = 'W1Re300He_111\peak745\depth';
    i= d745 + depth;
    s = strcat(str1,num2str(i));
    cd (s);
    Q3D = dlmread('Qspace3D.txt');
    QZ_coord = dlmread('Qspace3DCorners.txt');
    cd('..\..\..\');
    width = 4*(qwidth_745(pos_745+depth,2));
    if width ==0
        width = 1e-3;
    end
    power = 0.5;%8;%5
    qspace_Centre745_gauss(depth+1,:) = qspace_centre(Q3D,QZ_coord,width,power);
    
    %peak1024
    str1 = 'W1Re300He_111\peak1024\depth';
    i= d1024 + depth;
    s = strcat(str1,num2str(i));
    cd (s);
    Q3D = dlmread('Qspace3D.txt');
    QZ_coord = dlmread('Qspace3DCorners.txt');
    cd('..\..\..\');
    width = 4*(qwidth_1024(pos_1024+depth,2));
    if width ==0
        width = 1e-3;
    end
    power = 0.5;%8;%5
    qspace_Centre1024_gauss(depth+1,:) = qspace_centre(Q3D,QZ_coord,width,power);
    
    all_refl_qcenter (1:3,1) = qspace_Centre354_gauss(depth+1,1:3);
    all_refl_qcenter (1:3,2) = qspace_Centre444_gauss(depth+1,1:3);
    all_refl_qcenter (1:3,3) = qspace_Centre563_gauss(depth+1,1:3);
    all_refl_qcenter (1:3,4) = qspace_Centre804_gauss(depth+1,1:3);
    all_refl_qcenter (1:3,5) = qspace_Centre680_gauss(depth+1,1:3);
    all_refl_qcenter (1:3,6) = qspace_Centre745_gauss(depth+1,1:3);
    all_refl_qcenter (1:3,7) = qspace_Centre1024_gauss(depth+1,1:3);

    
    all_refl_hkl (1:3,1) = [3 5 4];
    all_refl_hkl (1:3,2) = [4 4 4];
    all_refl_hkl (1:3,3) = [5 6 3];
    all_refl_hkl (1:3,4) = [8 0 4];       
    all_refl_hkl (1:3,5) = [6 8 0];
    all_refl_hkl (1:3,6) = [7 4 5];
    all_refl_hkl (1:3,7) = [10 2 4];

    
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
%WITHOUT COMBINATION rotating by U
clear strain_depth_rot;

for i=1:29
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
    strain_depth_rot(i,i,1:27) = strain_depth_rot(i,i,1:27)+(0-mean(strain_depth_rot(i,i,18:27)));
%     strain_depth_rot(i,i,14:27)=0;
end

strain_depth_rot(1,2,1:27) = strain_depth_rot(1,2,1:27)+(0-mean(strain_depth_rot(1,2,23:27)));
strain_depth_rot(1,3,1:27) = strain_depth_rot(1,3,1:27)+(0-mean(strain_depth_rot(1,3,23:27)));
strain_depth_rot(2,3,1:27) = strain_depth_rot(2,3,1:27)+(0-mean(strain_depth_rot(2,3,23:27)));

%%
%the subtrate below is assumed to be strain free and thus being assigned
%the same value as at 8um depth
for i=28:31
    strain_depth_rot(:,:,i) = strain_depth_rot(:,:,27);
end

%%
%plotting principle strains without combination
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
clear strain;
strain(1,:) = squeeze(strain_depth_rot(1,2,1:31));
ylim([-lim,lim]);
xlim([0,13]);
p1 = plot(depth,strain);
p1.LineWidth = 2;
hold on;
xlabel('depth (\mum)','fontsize', 20,'FontWeight','bold');
ylabel('strain','fontsize', 20,'FontWeight','bold');
clear strain;
strain(1,:) = squeeze(strain_depth_rot(1,3,1:31));
p2 = plot(depth,strain);
ylim([-lim,lim]);
xlim([0,13]);
p2.LineWidth = 2;
hold on;
clear strain;
strain(1,:) = squeeze(strain_depth_rot(2,3,1:31));
p3 = plot(depth,strain);
ylim([-lim,lim]);
xlim([0,13]);
p3.LineWidth = 2;
h1 = legend('\epsilon_{xy}','\epsilon_{xz}','\epsilon_{yz}');
set(h1,'fontsize',13);
set(gca,'fontsize',18);
title(sprintf('Shear Strain'),'fontsize',15);



