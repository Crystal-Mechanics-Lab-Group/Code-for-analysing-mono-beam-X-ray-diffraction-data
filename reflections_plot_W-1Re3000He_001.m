%storing variables for given depth
B_depth = zeros(3,3,100);
U_depth = zeros(3,3,100);
T_depth = zeros(3,3,100);
A_depth = zeros(3,3,100);
strain_depth = zeros(3,3,100);
vol_strain_depth = zeros(3,3,100);
dev_strain_depth = zeros(3,3,100);

%surface of peaks at reconstructed depth as followss:
d47m3 = 12;
d181 = 12;
dm390 = 18;

cd('W1Re3000He_001');
%%

%reflections at depth
for depth=0:30

    %peak47-3
    str1 = 'peak47-3\depth';
    i= d47m3 + depth;
    s = strcat(str1,num2str(i));
    cd(s);
    Q3D = dlmread('Qspace3D.txt');
    QZ_coord = dlmread('Qspace3DCorners.txt');
    cd('..\..\');
    qspace_Centre47m3 = qspace_centre(Q3D,QZ_coord);
    
    %peak181
    str1 = 'peak181\depth';
    i= d181 + depth;
    s = strcat(str1,num2str(i));
    cd (s);
    Q3D = dlmread('Qspace3D.txt');
    QZ_coord = dlmread('Qspace3DCorners.txt');
    cd('..\..\');
    qspace_Centre181 = qspace_centre(Q3D,QZ_coord);
    
    %peak-390
    str1 = 'peak-390\depth';
    i= dm390 + depth;
    s = strcat(str1,num2str(i));
    cd (s);
    Q3D = dlmread('Qspace3D.txt');
    QZ_coord = dlmread('Qspace3DCorners.txt');
    cd('..\..\');
    qspace_Centrem390 = qspace_centre(Q3D,QZ_coord);
    
  
    all_refl_qcenter (1:3,1) = qspace_Centre47m3;
    all_refl_qcenter (1:3,2) = qspace_Centre181;
    all_refl_qcenter (1:3,3) = qspace_Centrem390;
    
    all_refl_hkl (1:3,1) = [4 7 -3];
    all_refl_hkl (1:3,2) = [1 8 1];
    all_refl_hkl (1:3,3) = [-3 9 0];

    
    %UB*all_refl_hkl = all_refl_qcenter OR UB =
    
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
        B_depth(:,:,i) = squeeze(B);
        U_depth(:,:,i) = squeeze(U);
        T_depth(:,:,i) = squeeze(T);
        A_depth(:,:,i) = squeeze(A);
        strain_depth(:,:,i) = squeeze(strain);
        vol_strain_depth(:,:,i) = squeeze(vol_strain);
        dev_strain_depth(:,:,i) = squeeze(dev_strain);
    end

end

%%

%rotate strain tensor using U_Crystal_to_Sample
for i=1:31
    rot = rotx(135)*U_depth(:,:,i);
    strain_depth_rot(:,:,i) = rot*strain_depth(:,:,i)*transpose(rot);
    vol_strain_depth_rot(:,:,i) = rot*vol_strain_depth(:,:,i)*transpose(rot);
    dev_strain_depth_rot(:,:,i) = rot*dev_strain_depth(:,:,i)*transpose(rot);
end

%%
%Considering the average of values below 7um depth to be 0, because it is
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

