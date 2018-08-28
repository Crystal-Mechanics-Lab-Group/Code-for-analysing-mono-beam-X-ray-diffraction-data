%storing variables for given depth
B_depth = zeros(3,3,100);
U_depth = zeros(3,3,100);
T_depth = zeros(3,3,100);
A_depth = zeros(3,3,100);
strain_depth = zeros(3,3,100);
vol_strain_depth = zeros(3,3,100);
dev_strain_depth = zeros(3,3,100);

%surface of peaks at reconstructed depth as followss:
d356 = 21;
d352 = 19;
dm372 = 15;

cd('W3000He_001');

%%

%reflections at depth
for depth=0:29

    %peak356
    str1 = 'peak356_2\depth';
    i= d356 + depth;
    s = strcat(str1,num2str(i));
    cd(s);
    Q3D = dlmread('Qspace3D.txt');
    QZ_coord = dlmread('Qspace3DCorners.txt');
    cd('..\..\');     
    qspace_Centre356(depth+1,:) = qspace_centre(Q3D,QZ_coord);
    
    %peak352
    str1 = 'peak352\depth';
    i= d352 + depth;
    s = strcat(str1,num2str(i));
    cd (s);
    Q3D = dlmread('Qspace3D.txt');
    QZ_coord = dlmread('Qspace3DCorners.txt');
    cd('..\..\');     
    qspace_Centre352(depth+1,:) = qspace_centre(Q3D,QZ_coord);
    
    %peak-372
    str1 = 'peak-372\depth';
    i= dm372 + depth;
    s = strcat(str1,num2str(i));
    cd (s);
    Q3D = dlmread('Qspace3D.txt');
    QZ_coord = dlmread('Qspace3DCorners.txt');
    cd('..\..\');     
    qspace_Centrem372(depth+1,:) = qspace_centre(Q3D,QZ_coord);
    

end


%%
for i=1:30
    all_refl_qcenter (1:3,1) = qspace_Centre356(i,:);
    all_refl_qcenter (1:3,2) = qspace_Centre352(i,:);
    all_refl_qcenter (1:3,3) = qspace_Centrem372(i,:);

    all_refl_hkl (1:3,1) = [3 5 6];
    all_refl_hkl (1:3,2) = [3 5 2];
    all_refl_hkl (1:3,3) = [-3 7 2];

    %UB*all_refl_hkl = all_refl_qcenter OR UB =

    UB = all_refl_qcenter/all_refl_hkl;
    [U, B] = UB2UandB_2(UB);
    [A, A0] = AfromB(B);

    %transformation matrix
    T = A/A0;
    I = [1 0 0; 0 1 0; 0 0 1];
    strain = 0.5*(T+transpose(T))-I;
    j = (1/3)*trace(strain);
    vol_strain = [j 0 0; 0 j 0; 0 0 j];
    dev_strain = strain - vol_strain;

    
    B_depth(:,:,i) = squeeze(B);
    U_depth(:,:,i) = squeeze(U);
    T_depth(:,:,i) = squeeze(T);
    A_depth(:,:,i) = squeeze(A);
    strain_depth(:,:,i) = squeeze(strain);
    vol_strain_depth(:,:,i) = squeeze(vol_strain);
    dev_strain_depth(:,:,i) = squeeze(dev_strain);
    
end

%%
%rotate strain tensor using U_Crystal_to_Sample
for i=1:30
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
    strain_depth_rot(i,i,1:30) = strain_depth_rot(i,i,1:30)+(0-mean(strain_depth_rot(i,i,23:30)));
%     strain_depth_rot(i,i,14:27)=0;
end

strain_depth_rot(1,2,1:30) = strain_depth_rot(1,2,1:30)+(0-mean(strain_depth_rot(1,2,23:30)));
strain_depth_rot(1,3,1:30) = strain_depth_rot(1,3,1:30)+(0-mean(strain_depth_rot(1,3,23:30)));
strain_depth_rot(2,3,1:30) = strain_depth_rot(2,3,1:30)+(0-mean(strain_depth_rot(2,3,23:30)));

%%
%plotting principle strains 
figure
depth = linspace(0,14.5,30);
depth = depth*cosd(45);
lim = 2e-3;
    for i=1:3
        clear strain;
        strain(1,:) = squeeze(strain_depth_rot(i,i,1:30));
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
depth = linspace(0,14.5,30);
depth = depth*cosd(45);
    clear strain;
    strain(1,:) = squeeze(strain_depth_rot(1,2,1:30));
    ylim([-lim,lim]);
    xlim([0,10]);
    p1 = plot(depth,strain);
    p1.LineWidth = 2;
    hold on;
    xlabel('depth (\mum)','fontsize', 20,'FontWeight','bold');
    ylabel('strain','fontsize', 20,'FontWeight','bold');
    clear strain;
    strain(1,:) = squeeze(strain_depth_rot(1,3,1:30));
    p2 = plot(depth,strain);
    ylim([-lim,lim]);
    xlim([0,10]);
    p2.LineWidth = 2;
    hold on;
    clear strain;
    strain(1,:) = squeeze(strain_depth_rot(2,3,1:30));
    p3 = plot(depth,strain);
    ylim([-lim,lim]);
    xlim([0,10]);
    p3.LineWidth = 2;
    h1 = legend('\epsilon_{xy}','\epsilon_{xz}','\epsilon_{yz}');
    set(h1,'fontsize',13);
    set(gca,'fontsize',18);
    title(sprintf('Shear Strain'),'fontsize',15);



