%storing variables for given combination at a given depth
B_depth = zeros(3,3,100,100);
U_depth = zeros(3,3,100,100);
T_depth = zeros(3,3,100,100);
A_depth = zeros(3,3,100,100);
strain_depth = zeros(3,3,100,100);
vol_strain_depth = zeros(3,3,100,100);
dev_strain_depth = zeros(3,3,100,100);

%surface of peaks at reconstructed depth as followss:
d039 = 10;
d008 = 7;
dm136 = 8;
d367 = 13;
d266 = 11;
dm217 = 5;


surf_039 = -5;
surf_008 = -6.5;
surf_m136 = -6;
surf_367 = -3.5;
surf_266 = -4.5;
surf_m217 = -7.5;

%%

cd('W300He_001');
qwidth_039 = dlmread('Q_Positions_Qwidth_039.txt');
qwidth_008 = dlmread('Q_Positions_Qwidth_008.txt');
qwidth_m136 = dlmread('Q_Positions_Qwidth_-136.txt');
qwidth_367 = dlmread('Q_Positions_Qwidth_367.txt');
qwidth_266 = dlmread('Q_Positions_Qwidth_266.txt');
qwidth_m217 = dlmread('Q_Positions_Qwidth_-217.txt');

%%
for i=1:size(qwidth_039,1)
    if qwidth_039(i,1) == surf_039
        pos_039 = i;
        break;
    end
end

for i=1:size(qwidth_008,1)
    if qwidth_008(i,1) == surf_008
        pos_008 = i;
        break;
    end
end

for i=1:size(qwidth_m136,1)
    if qwidth_m136(i,1) == surf_m136
        pos_m136 = i;
        break;
    end
end

for i=1:size(qwidth_367,1)
    if qwidth_367(i,1) == surf_367
        pos_367 = i;
        break;
    end
end

for i=1:size(qwidth_266,1)
    if qwidth_266(i,1) == surf_266
        pos_266 = i;
        break;
    end
end

for i=1:size(qwidth_m217,1)
    if qwidth_m217(i,1) == surf_m217
        pos_m217 = i;
        break;
    end
end

cd('..\');


%%
    
%reflections at depth
for depth=0:30

    %peak039
    str1 = 'W300He_001\peak039\depth';
    i= d039 + depth;
    s = strcat(str1,num2str(i));
    cd(s);
    Q3D = dlmread('Qspace3D.txt');
    QZ_coord = dlmread('Qspace3DCorners.txt');
    cd('..\..\..\');
    width = 4*(qwidth_039(pos_039+depth,2));
    if width ==0
        width = 1e-3;
    end
    power = 0.5;
    qspace_Centre039_gauss(depth+1,:) = qspace_centre(Q3D,QZ_coord,width,power);

   
   %peak008
    str1 = 'W300He_001\peak008\depth';
    i= d008 + depth;
    s = strcat(str1,num2str(i));
    cd(s);
    Q3D = dlmread('Qspace3D.txt');
    QZ_coord = dlmread('Qspace3DCorners.txt');
    cd('..\..\..\');
    width = 4*(qwidth_008(pos_008+depth,2));
    if width ==0
        width = 1e-3;
    end
    power = 0.5;
    qspace_Centre008_gauss(depth+1,:) = qspace_centre(Q3D,QZ_coord,width,power);


   %peakm136
    str1 = 'W300He_001\peak-136\depth';
    i= dm136 + depth;
    s = strcat(str1,num2str(i));
    cd(s);
    Q3D = dlmread('Qspace3D.txt');
    QZ_coord = dlmread('Qspace3DCorners.txt');
    cd('..\..\..\');
    width = 4*(qwidth_m136(pos_m136+depth,2));
    if width ==0
        width = 1e-3;
    end
    power = 0.5;
    qspace_Centrem136_gauss(depth+1,:) = qspace_centre(Q3D,QZ_coord,width,power);


   %peakm367
    str1 = 'W300He_001\peak367\depth';
    i= d367 + depth;
    s = strcat(str1,num2str(i));
    cd(s);
    Q3D = dlmread('Qspace3D.txt');
    QZ_coord = dlmread('Qspace3DCorners.txt');
    cd('..\..\..\');
    width = 4*(qwidth_367(pos_367+depth,2));
    if width ==0
        width = 1e-3;
    end
    power = 0.5;
    qspace_Centre367_gauss(depth+1,:) = qspace_centre(Q3D,QZ_coord,width,power);


    %peak266
    str1 = 'W300He_001\peak266\depth';
    i= d266 + depth;
    s = strcat(str1,num2str(i));
    cd(s);
    Q3D = dlmread('Qspace3D.txt');
    QZ_coord = dlmread('Qspace3DCorners.txt');
    cd('..\..\..\');
    width = 4*(qwidth_266(pos_266+depth,2));
    if width ==0
        width = 1e-3;
    end
    power = 0.5;
    qspace_Centre266_gauss(depth+1,:) = qspace_centre(Q3D,QZ_coord,width,power);

  
 %peak-217
    str1 = 'W300He_001\peak-217\depth';
    i= dm217 + depth;
    s = strcat(str1,num2str(i));
    cd(s);
    Q3D = dlmread('Qspace3D.txt');
    QZ_coord = dlmread('Qspace3DCorners.txt');
    cd('..\..\..\');
    width = 4*(qwidth_m217(pos_m217+depth,2));
    if width ==0
        width = 1e-3;
    end
    power = 0.5;
    qspace_Centrem217_gauss(depth+1,:) = qspace_centre(Q3D,QZ_coord,width,power);

    
      
    all_refl_qcenter (1:3,1) = qspace_Centre039_gauss(depth+1,1:3);
    all_refl_qcenter (1:3,2) = qspace_Centre008_gauss(depth+1,1:3);
    all_refl_qcenter (1:3,3) = qspace_Centrem136_gauss(depth+1,1:3);
    all_refl_qcenter (1:3,4) = qspace_Centre367_gauss(depth+1,1:3);
    all_refl_qcenter (1:3,5) = qspace_Centre266_gauss(depth+1,1:3);
    all_refl_qcenter (1:3,6) = qspace_Centrem217_gauss(depth+1,1:3);

    
    all_refl_hkl (1:3,1) = [0 3 9];
    all_refl_hkl (1:3,2) = [0 0 8];
    all_refl_hkl (1:3,3) = [-1 3 6];
    all_refl_hkl (1:3,4) = [3 6 7];       
    all_refl_hkl (1:3,5) = [2 6 6];
    all_refl_hkl (1:3,6) = [-2 1 7];

    
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
    strain_depth_rot(i,i,1:31) = strain_depth_rot(i,i,1:31)+(0-mean(strain_depth_rot(i,i,23:31)));
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
set(h1,'fontsize',18);
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



