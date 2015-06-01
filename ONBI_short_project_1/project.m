%% ONBI SHORT PROJECT 1
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Jack Allen
% Supervisor: Vicente Grau
%
% (DETERMINE = Myocardium infarction)
% (MESA = no symptoms of myocardium infarction)
% Initialise
clear all
close all
clc

setenv('path',[getenv('path'),';','C:\Users\jesu2687\Documents\MATLAB\ErnestoCode\Tools\MESHES\vtk_libs']);
addpath('C:\Users\jesu2687\Documents\MATLAB\ONBI_short_project_1')
addpath C:\Users\jesu2687\Documents\MATLAB\ErnestoCode\Tools\MESHES\
addpath C:\Users\jesu2687\Documents\MATLAB\ErnestoCode\Tools
addpath C:\Users\jesu2687\Documents\MATLAB\ErnestoCode\
addpath C:\Users\jesu2687\Documents\MATLAB\output-stacom-newcase\output-stacom-newcase
addpath C:\Users\jesu2687\Documents\MATLAB\closestPoint


load('C:\Users\jesu2687\Documents\MATLAB\ONBI_short_project_1\project1_data.mat');
load('C:\Users\jesu2687\Documents\MATLAB\ONBI_short_project_1\labels.mat');

%store indices that allow us to indentify which study (DETERMINE or MESA)
%each case belongs to.
data(1).DETERMINE_indices = Labels(2:101,3);
data(1).MESA_indices = Labels(102:201,3);
data(1).MESA_indices(2) = 401 ; % replace SMM0001 with the replacement provided by the organisers (SMM0401

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATE MYOCARDIUM THICKNESS (DISTANCE BETWEEN EPI AND ENDO)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('start myo thicknesses')
run('calc_myocardium_thicknesses.m')
disp('finished myo thicknesses')
%%
figure;
subplot 221
title 'DETERMINE - VARIANCE OF EACH DISTANCE (1089) OVER THE POPULATION'
patch( 'vertices',reshape(DETERMINEmeanDiaEpi,[1089,3]),'faces',data(401).diastolic.epi.tri,'facecolor','interp','cdata', 'edgecolor','none')%[1 1 1]*0.2)
axis equal;
view(3);
colormap jet
colorbar
caxis([-3 7])
subplot 222
title 'MESA - '
patch( 'vertices',reshape(DETERMINEmeanDiaEpi,[1089,3]),'faces',data(401).diastolic.epi.tri,'facecolor','interp','cdata',  'edgecolor','none')%[1 1 1]*0.2)
axis equal;
view(3);
colormap jet
colorbar
caxis([-3 7])


figure;
subplot 221
title 'DETERMINE - mean myo thickness change'
patch( 'vertices',reshape(DETERMINEmeanDiaEpi,[1089,3]),'faces',data(401).diastolic.epi.tri,'facecolor','interp','cdata',mean(cell2mat({data(data(1).DETERMINE_indices).myo_T_changes}),2),'edgecolor','none')%[1 1 1]*0.2)
axis equal;
view(3);
colormap jet
colorbar
caxis([-3 7])
subplot 222
title 'MESA - mean myo thickness change'
patch( 'vertices',reshape(DETERMINEmeanDiaEpi,[1089,3]),'faces',data(401).diastolic.epi.tri,'facecolor','interp','cdata',mean(cell2mat({data(data(1).MESA_indices).myo_T_changes}),2),'edgecolor','none')%[1 1 1]*0.2)
axis equal;
view(3);
colormap jet
colorbar
caxis([-3 7])


figure;
subplot 221
title 'DETERMINE - MEAN DIASTOLE'
patch( 'vertices',reshape(DETERMINEmeanDiaEpi,[1089,3]),'faces',data(401).diastolic.epi.tri,'facecolor','interp','cdata',DETERMINEmeanDiaT,'edgecolor','none')%[1 1 1]*0.2)
axis equal;
view(3);
colormap jet
colorbar
caxis([1 16])
subplot 222
title 'MESA - MEAN DIASTOLE'
patch( 'vertices',reshape(MESAmeanDiaEpi,[1089,3]),'faces',data(401).diastolic.epi.tri,'facecolor','interp','cdata',MESAmeanDiaT,'edgecolor','none')%[1 1 1]*0.2)
axis equal;
view(3);
colormap jet
colorbar
caxis([1 16])

subplot 223
title 'DETERMINE - MEAN SYSTOLE'
patch( 'vertices',reshape(DETERMINEmeanSysEpi,[1089,3]),'faces',data(401).systolic.epi.tri,'facecolor','interp','cdata',DETERMINEmeanSysT,'edgecolor','none')%[1 1 1]*0.2)
axis equal;
view(3);
colormap jet
colorbar
caxis([1 16])
subplot 224
title 'MESA - MEAN SYSTOLE'
patch( 'vertices',reshape(MESAmeanSysEpi,[1089,3]),'faces',data(401).systolic.epi.tri,'facecolor','interp','cdata',MESAmeanSysT,'edgecolor','none')%[1 1 1]*0.2)
axis equal;
view(3);
colormap jet
colorbar
caxis([1 16])
%% plot myo thickness stats
figure
subplot 121
hold on
histogram(cell2mat({deltaMeanT(data(1).MESA_indices)}))
histogram(cell2mat({deltaMeanT(data(1).DETERMINE_indices)}))
xlabel 'change in mean thickness from systole to diastole'
hold off
[data, accuracy, sensitivity, specificity] = calcAccuracy( data, cell2mat({deltaMeanT(data(1).DETERMINE_indices)}),  cell2mat({deltaMeanT(data(1).MESA_indices)}), 1,8,1);
subplot 122
hold on
plot(accuracy)
plot(sensitivity)
plot(specificity)

figure
subplot 231
title 'diastolic, mean myocardium thickness'
hold on
histogram(cell2mat({diaMeanT(data(1).MESA_indices)}))
histogram(cell2mat({diaMeanT(data(1).DETERMINE_indices)}))
ylabel 'number of LVs'
xlabel 'mean'
legend 'MESA' 'DETERMINE'
subplot 232
title 'systolic, mean myocardium thickness'
hold on
histogram(cell2mat({sysMeanT(data(1).MESA_indices)}))
histogram(cell2mat({sysModeT(data(1).DETERMINE_indices)}))
subplot 233
title 'diastolic, mode myocardium thickness'
hold on
histogram(cell2mat({diaModeT(data(1).MESA_indices)}))
histogram(cell2mat({diaModeT(data(1).DETERMINE_indices)}))
subplot 234
title 'systolic, mode myocardium thickness'
hold on
histogram(cell2mat({sysModeT(data(1).MESA_indices)}))
histogram(cell2mat({sysModeT(data(1).DETERMINE_indices)}))
subplot 235
title 'diastolic, median myocardium thickness'
hold on
histogram(cell2mat({diaMedianT(data(1).MESA_indices)}))
histogram(cell2mat({diaMedianT(data(1).DETERMINE_indices)}))
subplot 236
title 'systolic, median myocardium thickness'
hold on
histogram(cell2mat({sysMedianT(data(1).MESA_indices)}))
histogram(cell2mat({sysMedianT(data(1).DETERMINE_indices)}))

figure
subplot 121
hold on
histogram(cell2mat({data(data(1).MESA_indices).dia_dEPI2ENDO}))
histogram(cell2mat({data(data(1).DETERMINE_indices).dia_dEPI2ENDO}))
legend 'MESA' 'DETERMINE'
xlabel 'dEPI2ENDO'
title 'diastolic myocardium thicknesses'
subplot 122
hold on
histogram(cell2mat({data(data(1).MESA_indices).sys_dEPI2ENDO}))
histogram(cell2mat({data(data(1).DETERMINE_indices).sys_dEPI2ENDO}))
legend 'MESA' 'DETERMINE'
xlabel 'dEPI2ENDO'
title 'systolic myocardium thicknesses'

figure
hold on
histogram(cell2mat({data(data(1).MESA_indices).myo_T_changes}))
histogram(cell2mat({data(data(1).DETERMINE_indices).myo_T_changes}))
legend 'MESA' 'DETERMINE'
xlabel 'myocardium thickness change from systole to diastole'

figure
subplot 221
hold on
histogram(cell2mat({data(data(1).MESA_indices).dia_dEPI2ENDO_vars}))
histogram(cell2mat({data(data(1).DETERMINE_indices).dia_dEPI2ENDO_vars}))
legend 'MESA' 'DETERMINE'
xlabel 'dia dEPI2ENDO variances'
subplot 222
hold on
histogram(cell2mat({data(data(1).MESA_indices).sys_dEPI2ENDO_vars}))
histogram(cell2mat({data(data(1).DETERMINE_indices).sys_dEPI2ENDO_vars}))
legend 'MESA' 'DETERMINE'
xlabel 'sys dEPI2ENDO variances'
subplot 223
hold on
histogram(cell2mat({data(data(1).MESA_indices).dia_dEPI2ENDO_stds}))
histogram(cell2mat({data(data(1).DETERMINE_indices).dia_dEPI2ENDO_stds}))
legend 'MESA' 'DETERMINE'
xlabel 'dia dEPI2ENDO standard deviations'
subplot 224
hold on
histogram(cell2mat({data(data(1).MESA_indices).sys_dEPI2ENDO_stds}))
histogram(cell2mat({data(data(1).DETERMINE_indices).sys_dEPI2ENDO_stds}))
legend 'MESA' 'DETERMINE'
xlabel 'sys dEPI2ENDO standard deviations'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATE MYOCARDIUM VOLUMES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% when run after calcVolumes() etc,this gives warnings...

disp('start making myo shapes')

%load struct containing the triangle file for the myo 'donut' shaped lid
load('C:\Users\jesu2687\Documents\MATLAB\ONBI_short_project_1\myoB_tri.mat')

% make myo shapes (using code from 'calcMyoVolume()' )
  [data] = calcMyoVolume(data, myoB);

disp('finished making myo shapes')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATE ENDO AND EPI VOLUMES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('calculating endo and epi volumes')

%calculate volumes
[data, sys_epi_volumes, sys_endo_volumes, dia_epi_volumes, dia_endo_volumes] = calcVolumes(data);
%store volumes
[data] = storeVolumes(data, sys_epi_volumes, sys_endo_volumes, dia_epi_volumes, dia_endo_volumes);

disp('finished calculating endo and epi volumes')

plotVolumes(data,sys_epi_volumes,sys_endo_volumes,dia_epi_volumes,dia_endo_volumes)

%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATE MYO VOLUMES (EPI - ENDO)
disp('start EPI - ENDO')
[data] = epiMinusEndoVolumes(data);
disp('finish EPI - ENDO')
%% %%%%%%%%%%%%%%%%%%%%%%%

% store volumes for plotting
[data] = storeMyoVolumes(data);

%% compare my myo volumes with Ernesto's
ErnestoVolumes = importdata('volumes.txt');
figure
xlim ([0 10e5]);
ylim = ([0 10e5]);
hold on
plot(cell2mat({data([data(1).DETERMINE_indices ; data(1).MESA_indices]).dia_myovolumes}), ErnestoVolumes.data(:,4), '.');
plot(cell2mat({data([data(1).DETERMINE_indices ; data(1).MESA_indices]).sys_myovolumes}), ErnestoVolumes.data(:,3), '.');
title 'myocardium volumes'
xlabel 'Jack volumes'
ylabel 'Ernesto volumes'
legend 'diastolic' 'systolic'

%import volumes that Ernesto calculated
ErnestoVolumes = importdata('volumes.txt');
figure
xlim ([0 10e5]);
ylim = ([0 10e5]);
hold on
plot(cell2mat({data([data(1).DETERMINE_indices ; data(1).MESA_indices]).dia_epiMinusEndoVolumes}), ErnestoVolumes.data(:,4), '.');
plot(cell2mat({data([data(1).DETERMINE_indices ; data(1).MESA_indices]).sys_epiMinusEndoVolumes}), ErnestoVolumes.data(:,3), '.');
title 'myocardium volumes'
xlabel 'Jack volumes'
ylabel 'Ernesto volumes'
legend 'diastolic' 'systolic'

% plot systolic myocardium volumes
plotMyoVolumes(data)

%% calculate "myocardium ejection fraction" myoEF
% myoSV = diastolic myovolume - systolic myovolume
% myoEF = (myoSV/(diastolic volume))*100
[data] = calcMyoEjectionFractions(data);
%plot myocardium "ejection fractions"
figure
nBins = 25;
title 'myocardium "ejection fractions"'
hold on
histogram(cell2mat({data((data(1).DETERMINE_indices)').myoEF}), nBins)
histogram(cell2mat({data((data(1).MESA_indices)').myoEF}), nBins)
xlabel 'Ejection fraction (%)'
ylabel Frequency

%%
% use Ernesto's myo volumes
% for i = 1:401
%
%     data(i).

% for i = 1:401
% all(i,:) = data(i).dia_dEPI2ENDO;
% end
%    mean = mean(all);
%
% figure; patch( 'vertices',data(1).diastolic.endo.xyz,'faces',data(1).diastolic.epi.tri,'facecolor','interp','cdata',mean','edgecolor',[1 1 1]*0.2)
% axis equal;
% view(3);
% colormap jet
% colorbar

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  CALCULATE EJECTION FRACTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('started calculating ejection fractions')

[data, DETERMINE_strokeVolumes, DETERMINE_ejectionFractions,  MESA_strokeVolumes,  MESA_ejectionFractions] = calcEjectionFraction(data);

plotEF(data)
disp('finished calculating ejection fractions')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  CALCULATE ACCURACIES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('started calculating accuracies')

%pass the positive and negative condition
[data, accuracyEF, sensitivityEF, specificityEF] = calcAccuracy( data, cell2mat({data(data(1).DETERMINE_indices).ejectionFraction}) ,  cell2mat({data(data(1).MESA_indices).ejectionFraction}), 1,100,1);

figure
plotROC(sensitivityEF, specificityEF)
title 'Ejection fraction'

plotAccuracyEF(data, accuracyEF, cell2mat({data(data(1).DETERMINE_indices).ejectionFraction}) ,  cell2mat({data(data(1).MESA_indices).ejectionFraction}), 'DETERMINE', 'MESA' , 100)
xlabel 'ejection fraction (%)'
ylabel 'frequency'

[data, data(1).accuracies.EF, data(1).sensitivities.EF, data(1).specificities.EF] = calcThresholdAccuracy(data, cell2mat({data(data(1).DETERMINE_indices).ejectionFraction}) ,  cell2mat({data(data(1).MESA_indices).ejectionFraction}),1,54)

disp('finished calculating accuracies')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATE MESH AREAS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gives same numbers as Ernesto's code :) (MeshAreatriangles(M))
disp('started calculating triangular mesh areas and triangle side lengths')
[data] = calcTriMeshAreas(data);

disp('finished calculating triangular mesh areas and triangle side lengths')

% calculate area to volume ratios
[data] = calcAVR(data);

DETERMINE_sys_endo_AVratios = cell2mat({data(:).DETERMINE_sys_endo_AVratio})' ;
DETERMINE_dia_endo_AVratios = cell2mat({data(:).DETERMINE_dia_endo_AVratio})' ;
MESA_sys_endo_AVratios = cell2mat({data(:).MESA_sys_endo_AVratio})' ;
MESA_dia_endo_AVratios = cell2mat({data(:).MESA_dia_endo_AVratio})' ;

figure
subplot 221
nbins = 30;
title 'systolic myocardium, surface area to volume ratio'
hold on
histogram(cell2mat({data(data(1).MESA_indices).MESA_sys_myo_AVratio}),nbins)
histogram(cell2mat({data(data(1).DETERMINE_indices).DETERMINE_sys_myo_AVratio}),nbins)
xlabel 'surface area to volume ratio (AVR)'
ylabel 'frequency'
subplot 222
nbins = 30;
title 'systolic endocardium, surface area to volume ratio'
hold on
histogram(cell2mat({data(data(1).MESA_indices).MESA_sys_endo_AVratio}),nbins)
histogram(cell2mat({data(data(1).DETERMINE_indices).DETERMINE_sys_endo_AVratio}),nbins)
xlabel 'endocardium surface area to volume ratio (AVR)'
ylabel 'frequency'
subplot 223
nbins = 30;
title 'systolic epicardium, surface area to volume ratio'
hold on
histogram(cell2mat({data(data(1).MESA_indices).MESA_sys_epi_AVratio}),nbins)
histogram(cell2mat({data(data(1).DETERMINE_indices).DETERMINE_sys_epi_AVratio}),nbins)
xlabel 'epicardium surface area to volume ratio (AVR)'
ylabel 'frequency'
legend ' MESA' ' DETERMINE'

figure
subplot 221
nbins = 20;
title 'diastolic myocardium, surface area to volume ratio'
hold on
histogram(cell2mat({data(data(1).MESA_indices).MESA_dia_myo_AVratio}),nbins)
histogram(cell2mat({data(data(1).DETERMINE_indices).DETERMINE_dia_myo_AVratio}),nbins)
xlabel 'surface area to volume ratio (AVR)'
ylabel 'frequency'
subplot 222
nbins = 30;
title 'diastolic endocardium, surface area to volume ratio'
hold on
histogram(cell2mat({data(data(1).MESA_indices).MESA_dia_endo_AVratio}),nbins)
histogram(cell2mat({data(data(1).DETERMINE_indices).DETERMINE_dia_endo_AVratio}),nbins)
xlabel 'endocardium surface area to volume ratio (AVR)'
ylabel 'frequency'
subplot 223
nbins = 30;
title 'diastolic epicardium, surface area to volume ratio'
hold on
histogram(cell2mat({data(data(1).MESA_indices).MESA_dia_epi_AVratio}),nbins)
histogram(cell2mat({data(data(1).DETERMINE_indices).DETERMINE_dia_epi_AVratio}),nbins)
xlabel 'epicardium surface area to volume ratio (AVR)'
ylabel 'frequency'
legend ' MESA' 'DETERMINE'

[data, accuracyAV, sensitivityAV, specificityAV] = calcAccuracy(data, cell2mat({data(data(1).DETERMINE_indices).DETERMINE_sys_endo_AVratio}), cell2mat({data(data(1).MESA_indices).MESA_sys_endo_AVratio}),500,700,100);
figure
plotROC(sensitivityAV, specificityAV)
title 'surface area to volume ratio'
% legend ('systolic, endocardium surface area to volume ratio' , 'ejection fraction' ,'Location', 'best')

plotAccuracyAV( data,accuracyAV, cell2mat({data(data(1).DETERMINE_indices).DETERMINE_sys_endo_AVratio}), cell2mat({data(data(1).MESA_indices).MESA_sys_endo_AVratio}), 'DETERMINE', 'MESA', 25)
subplot 121
xlabel 'systolic endocardium surface area to volume ratio (%)'
subplot 122
xlabel 'systolic endocardium surface area to volume ratio (%)'

[data, data(1).accuracies.sysEndoAVratio, data(1).sensitivities.sysEndoAVratio, data(1).specificities.sysEndoAVratio] = calcThresholdAccuracy(data,cell2mat({data(data(1).DETERMINE_indices).DETERMINE_sys_endo_AVratio}) ,  cell2mat({data(data(1).MESA_indices).MESA_sys_endo_AVratio}),100,15)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GENERAL PROCRUSTES ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is done by concatenating the endo and epi shape vectors to produce
% longer vectors (myo.xyz). Procrustes analysis is then performed on the
% two myo.xyz vectors. Finally, the transformed myo.xyz vectors are split
% to extract transformed endo and epi shape vectors.
disp('starting procrustes analysis')
tic

%store dimensions of the shapes
[shape_nRows , shape_nCols] = size(data(2).diastolic.endo.xyz(1:1089,:));

%use the initial data clouds as references
initial_reference_id = 1; %sets which case will be the initial reference

%remove the B points (added in calcVolumes())
dia_endo_reference = data(initial_reference_id).diastolic.endo.xyz(1:1089,:);
dia_epi_reference = data(initial_reference_id).diastolic.epi.xyz(1:1089,:);
sys_endo_reference = data(initial_reference_id).systolic.endo.xyz(1:1089,:);
sys_epi_reference = data(initial_reference_id).systolic.epi.xyz(1:1089,:);

dia_myo_reference = reshape([dia_endo_reference(:) ; dia_epi_reference(:)], [2*1089 3]);
sys_myo_reference = reshape([sys_endo_reference(:)  ; sys_epi_reference(:)], [2*1089 3]);

dia_sys_myo_reference = [dia_myo_reference ; sys_myo_reference ];

% Think of endo and epi as one shape (concatenate).
for i = 1:401 %SMM001 has already been replaced by SMM401 in the script
    data(i).diastolic.myo.xyz = [data(i).diastolic.endo.xyz(1:1089,:) ; data(i).diastolic.epi.xyz(1:1089,:)];
    data(i).systolic.myo.xyz = [data(i).systolic.endo.xyz(1:1089,:) ; data(i).systolic.epi.xyz(1:1089,:)];
    
    %dia endo ; dia epi ; sys endo ; sys epi
    data(i).dia_sys.myo.xyz = [data(i).diastolic.myo.xyz ; data(i).systolic.myo.xyz] ;
end

%p = number of times procrustes is performed.
procrustes_iterations = 3;
[data, dia_myo_mean, sys_myo_mean, dia_sys_myo_mean, all_training_diastolic_myo_shapes, all_training_systolic_myo_shapes , all_training_dia_sys_myo_shapes] = calcProcrustes(data, procrustes_iterations, dia_myo_reference, sys_myo_reference, dia_sys_myo_reference);

% % reshape individual endo and epi shapes.
% dia_endo_mean = reshape(dia_endo_mean,size(data(1).diastolic.endo.xyz));
% dia_epi_mean = reshape(dia_epi_mean,size(data(1).diastolic.endo.xyz));
% sys_endo_mean = reshape(sys_endo_mean,size(data(1).diastolic.endo.xyz));
% sys_epi_mean = reshape(sys_epi_mean,size(data(1).diastolic.endo.xyz));

% dia_myo_reference = reshape(dia_myo_mean, [2*shape_nRows, shape_nCols]);
% sys_myo_reference = reshape(sys_myo_mean, [2*shape_nRows, shape_nCols]);

[data] = extractEndoEpi(data, shape_nRows);

disp('finished procrustes analysis')
toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SHAPE MODELING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PCA
% see Cootes tutorial for help: http://personalpages.manchester.ac.uk/staff/timothy.f.cootes/Models/app_models.pdf
disp('started PCA')
tic
% nCases = 400;
% for i = 1:nCases

% do PCA on MESA and DETERMINE cases only
for d = sort(data(1).DETERMINE_indices')
    
    DETERMINE_diastolic_myo_shapes(d,:) = all_training_diastolic_myo_shapes(d,:);
    DETERMINE_systolic_myo_shapes(d,:) = all_training_systolic_myo_shapes(d,:);
    
    DETERMINE_dia_sys_myo_shapes(d,:) = all_training_dia_sys_myo_shapes(d,:);
end
for m = sort(data(1).MESA_indices')
    MESA_diastolic_myo_shapes(m,:) = all_training_diastolic_myo_shapes(m,:);
    MESA_systolic_myo_shapes(m,:) = all_training_systolic_myo_shapes(m,:);
    
    MESA_dia_sys_myo_shapes(m,:) = all_training_dia_sys_myo_shapes(m,:);
end
% delete empty rows
DETERMINE_diastolic_myo_shapes( ~any(DETERMINE_diastolic_myo_shapes,2), : ) = [];
DETERMINE_systolic_myo_shapes( ~any(DETERMINE_systolic_myo_shapes,2), : ) = [];
DETERMINE_dia_sys_myo_shapes( ~any(DETERMINE_dia_sys_myo_shapes,2), : ) = [];
MESA_diastolic_myo_shapes( ~any(MESA_diastolic_myo_shapes,2), : ) = [];
MESA_systolic_myo_shapes( ~any(MESA_systolic_myo_shapes,2), : ) = [];
MESA_dia_sys_myo_shapes( ~any(MESA_dia_sys_myo_shapes,2), : ) = [];

%% check procrustes
figure
for m = 1:100
    m
    subplot 121
    title 'MESA, all four shapes'
    plot3D( reshape(MESA_dia_sys_myo_shapes(m,:),[4*1089 3]) ,'o')
    hold on
    
    subplot 122
    title 'DETERMINE, all four shapes'
    plot3D( reshape(DETERMINE_dia_sys_myo_shapes(m,:),[4*1089 3]) ,'o')
    hold on
    
end
%%  set shapes to train classifier
training_diastolic_myo_shapes = [DETERMINE_diastolic_myo_shapes ; MESA_diastolic_myo_shapes];
training_systolic_myo_shapes = [DETERMINE_systolic_myo_shapes ; MESA_systolic_myo_shapes];
training_dia_sys_myo_shapes = [DETERMINE_dia_sys_myo_shapes ; MESA_dia_sys_myo_shapes];

disp('started calculating covariance and mean shape')

dia_myo_cov_mat = cov(training_diastolic_myo_shapes);
sys_myo_cov_mat = cov(training_systolic_myo_shapes);
dia_sys_myo_cov_mat = cov(training_dia_sys_myo_shapes);

disp('finished calculating covariance and mean shape')
%% PCA - find eigenvectors of covariance matrix (of just MESA)
disp('started calculating eigenvectors of covariance matrix')
tic
[dia_myo_eigenvectors,dia_myo_eigenvalues]=eig(dia_myo_cov_mat);
[sys_myo_eigenvectors,sys_myo_eigenvalues]=eig(sys_myo_cov_mat);
[dia_sys_myo_eigenvectors,dia_sys_myo_eigenvalues]=eig(dia_sys_myo_cov_mat);
toc
% % Each eigenvalue gives the variance of the data about the mean in the
% % direction of the corresponding eigenvector. Compute the total variance
%  totalVariance = sum(dia_endo_eigenvalues(:));
% proportion = 1;
% % %Choose the first t largest eigenvalues such that
%  p = proportion*totalVariance;
% a = dia_endo_eigenvalues/sum(dia_endo_eigenvalues(:));
disp('finished calculating eigenvectors of covariance matrix')
%
nEigenvalues = 10;
sorted_dia_myo_eigenvalues = sort(dia_myo_eigenvalues(:), 'descend');
sorted_sys_myo_eigenvalues = sort(sys_myo_eigenvalues(:), 'descend');
sorted_dia_sys_myo_eigenvalues = sort(dia_sys_myo_eigenvalues(:), 'descend');

principle_dia_myo_eigenvalues = sorted_dia_myo_eigenvalues(1:nEigenvalues);
principle_sys_myo_eigenvalues = sorted_sys_myo_eigenvalues(1:nEigenvalues);
principle_dia_sys_myo_eigenvalues = sorted_dia_sys_myo_eigenvalues(1:nEigenvalues);

sorted_dia_myo_eigenvalues_contributions = 100*(sorted_dia_myo_eigenvalues./(sum(sorted_dia_myo_eigenvalues)));
sorted_sys_myo_eigenvalues_contributions = 100*(sorted_sys_myo_eigenvalues./(sum(sorted_sys_myo_eigenvalues)));
sorted_dia_sys_myo_eigenvalues_contributions = 100*(sorted_dia_sys_myo_eigenvalues./(sum(sorted_dia_sys_myo_eigenvalues)));

%% PCA - contribution from the chosen eigenvectors, as a percentage
principle_dia_myo_eigenvalues_contribution = sorted_dia_myo_eigenvalues_contributions(1:nEigenvalues);
principle_sys_myo_eigenvalues_contribution = sorted_sys_myo_eigenvalues_contributions(1:nEigenvalues);
principle_dia_sys_myo_eigenvalues_contribution = sorted_dia_sys_myo_eigenvalues_contributions(1:nEigenvalues);

principle_dia_sys_myo_eigenvectors = zeros(13068, nEigenvalues);
principle_dia_myo_eigenvectors = zeros(6534, nEigenvalues);
principle_sys_myo_eigenvectors = zeros(6534, nEigenvalues);

for n = 1:nEigenvalues
    [dia_myo_eRows(1,n), dia_myo_eCols(1,n)] = find(dia_myo_eigenvalues == principle_dia_myo_eigenvalues(n,1));
    [sys_myo_eRows(1,n), sys_myo_eCols(1,n)] = find(sys_myo_eigenvalues == principle_sys_myo_eigenvalues(n,1));
    [dia_sys_myo_eRows(1,n), dia_sys_myo_eCols(1,n)] = find(dia_sys_myo_eigenvalues == principle_dia_sys_myo_eigenvalues(n,1));
    
    principle_dia_myo_eigenvectors(:,n) = dia_myo_eigenvectors(:, dia_myo_eCols(1,n));
    principle_sys_myo_eigenvectors(:,n) = sys_myo_eigenvectors(:, sys_myo_eCols(1,n));
    principle_dia_sys_myo_eigenvectors(:,n) = dia_sys_myo_eigenvectors(:, dia_sys_myo_eCols(1,n)); 
end

%% Find model parameter values (b) , using +/- sqrt(eigenvalue) as b
% b can be increased later.

% systolic, endocardium
dia_myo_min_b = - sqrt(principle_dia_myo_eigenvalues);
dia_myo_max_b = sqrt(principle_dia_myo_eigenvalues);
sys_myo_min_b = -sqrt(principle_sys_myo_eigenvalues);
sys_myo_max_b = sqrt(principle_sys_myo_eigenvalues);
dia_sys_myo_min_b = - sqrt(principle_dia_sys_myo_eigenvalues);
dia_sys_myo_max_b = sqrt(principle_dia_sys_myo_eigenvalues);

disp('finished PCA')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Point distribution model (PDM)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('started point distribution model')
%% PDM - find mean shapes
mean_diastolic_myo_shape = mean(training_diastolic_myo_shapes);
mean_systolic_myo_shape = mean(training_systolic_myo_shapes);
mean_dia_sys_myo_shape = mean(training_dia_sys_myo_shapes);

%% PDM - find b values for the DETERMINE and MESA cases
for b = 1:10
    for i = data(1).DETERMINE_indices'
        for d = find(data(1).DETERMINE_indices==i)
            
            DETERMINE_dia_myo_b(d,b) = (principle_dia_myo_eigenvectors(:,b)')*(data(i).diastolic.myo.xyz(:) - mean_diastolic_myo_shape');
            DETERMINE_sys_myo_b(d,b) = (principle_sys_myo_eigenvectors(:,b)')*(data(i).systolic.myo.xyz(:) - mean_systolic_myo_shape');
            DETERMINE_dia_sys_myo_b(d,b) = (principle_dia_sys_myo_eigenvectors(:,b)')*(data(i).dia_sys.myo.xyz(:) - mean_dia_sys_myo_shape');
            
        end
    end
    
    for i = data(1).MESA_indices'
        for d = find(data(1).MESA_indices==i)
            %         d
            MESA_dia_myo_b(d,b) = (principle_dia_myo_eigenvectors(:,b)')*(data(i).diastolic.myo.xyz(:) - mean_diastolic_myo_shape');
            MESA_sys_myo_b(d,b) = (principle_sys_myo_eigenvectors(:,b)')*(data(i).systolic.myo.xyz(:) - mean_systolic_myo_shape');
            MESA_dia_sys_myo_b(d,b) = (principle_dia_sys_myo_eigenvectors(:,b)')*(data(i).dia_sys.myo.xyz(:) - mean_dia_sys_myo_shape');
            
        end
    end
end
disp('finished point distribution model')
%% visualise modes

%% visualise eigenvalue contributions
dia_myo_normalisedEigvalues = sorted_dia_myo_eigenvalues/(sum(sorted_dia_myo_eigenvalues));
sys_myo_normalisedEigvalues = sorted_sys_myo_eigenvalues/(sum(sorted_sys_myo_eigenvalues));
dia_sys_myo_normalisedEigvalues = sorted_dia_sys_myo_eigenvalues/(sum(sorted_dia_sys_myo_eigenvalues));
dia_sys_CS = cumsum(dia_sys_myo_normalisedEigvalues);
dia_CS = cumsum(dia_myo_normalisedEigvalues);
sys_CS = cumsum(sys_myo_normalisedEigvalues);

figure

hold on
plot(dia_CS(1:100))
plot(sys_CS(1:100))
for eigVal = 1:10
    plot(eigVal,(dia_CS(eigVal)),'o')
end
for eigVal = 1:10
    plot(eigVal,(sys_CS(eigVal)),'*')
end
xlabel 'eigenmode'
ylabel 'contribution'
legend 'diastole' 'systole'
hold off

figure
hold on
plot(dia_sys_CS(1:100))
for eigVal = 1:10
    plot(eigVal,(dia_sys_CS(eigVal)),'o')
end
xlabel 'eigenmode'
ylabel 'contribution'
legend 'diastole and systole'
hold off

%% visualise mode variations
visualisingMode = 2
stage = 'sys'
%animation
for n= 1
    
    x = [-60 ; 50];
    y = [-55 ; 50];
    z = [-110 ; 15];
    loops = 10;
    F(loops) = struct('cdata',[],'colormap',[]);
    
    for f = 1:10
        tmpf = sort(1:10,'descend')
        c = tmpf(f); % -1, -0.9, ... , -0.1
        pause(0.1)
        
        visualiseModesMovie(data, dia_sys_myo_mean,principle_dia_sys_myo_eigenvectors, dia_sys_myo_max_b, stage , visualisingMode, c,x,y,z)
        drawnow
        F(f) = getframe;
    end
    for f = 1:10
        
        c = -f % 0.1, 0.2, ... , 1
        pause(0.1)
        hold off
        visualiseModesMovie(data, dia_sys_myo_mean,principle_dia_sys_myo_eigenvectors, dia_sys_myo_max_b, stage, visualisingMode, c,x,y,z)
        drawnow
        F(f+10) = getframe;
    end
end

%% visualise modes of variation - static
% close all
visualisingMode1 = 1
visualisingMode2 = 2
c = 3 %number of standard deviations
x = [-60 ; 40];
y = [-65 ;50];
z = [-110 ; 20];
figure('name','diastolic modes of variation')
run('visualisemodes_diastole.m')
run('visualisemodes_diastole__combined_modes.m')
figure('name','systolic modes of variation')
run('visualisemodes_systole.m')
run('visualisemodes_systole__combined_modes.m')

% visualiseModes(data, dia_sys_myo_mean, principle_dia_sys_myo_eigenvectors, dia_sys_myo_max_b, 'sys', eMode, c,x,y,z)
%% play mode of variation movie
figure
title ([ 'mode',  num2str(visualisingMode) ])
axis off
for n = 1:4
    movie(F)
    
end
%% visualise b values - histograms
b = 2;

for b = 1:10
    b
    figure
    nbins = 30;
    hold on
    histogram(DETERMINE_sys_myo_b(:,b), nbins)
    histogram(MESA_sys_myo_b(:,b), nbins)
    hold off
    pause
end

% visualise b values - 2D plot
b = [ 2;5 ] ;
figure
hold on
xlabel 'b1'
ylabel 'b2'
plot(DETERMINE_sys_myo_b(:,b(1)), DETERMINE_sys_myo_b(:,b(2)), 'o')
plot(MESA_sys_myo_b(:,b(1)), MESA_sys_myo_b(:,b(2)), 'o')
hold off

% visualise b values - 3D plot
b = 1:3;
figure
hold on
plot3D(DETERMINE_sys_myo_b(:,b))
plot3D(MESA_sys_myo_b(:,b))
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Bullseye plots

figure('name','myocardium thickness bullseye plot: mode of variation , -/+*c*std')
colormap jet
cmin = 0;
cmax = 15;
subplot 121
axis tight
axis square
tmpShape = reshape(dia_sys_myo_new_shape_minus(:,visualisingMode1), [2*2178 , 3]);
shape = tmpShape(1:2178,:);
run('preProcessData_BullsEyePlots_NewShape.m');
subplot 122
axis tight
axis square
tmpShape = reshape(dia_sys_myo_new_shape_plus(:,visualisingMode1), [2*2178 , 3]);
shape = tmpShape(1:2178,:);
run('preProcessData_BullsEyePlots_NewShape.m')

figure('Name','Myocardium thickness bullseye plot - Mean') %,'NumberTitle','off')
colormap jet
axis tight
axis square
shape = reshape(dia_sys_myo_mean(:), [2*2178 , 3]);
shape = shape(1:2178,:);
run('preProcessData_BullsEyePlots_NewShape.m');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CLASSIFICATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SVM - set training data
trainingdata(:,1) = [DETERMINE_ejectionFractions' ; MESA_ejectionFractions'];
%modes of variation
trainingdata(:,2:6) = [DETERMINE_sys_myo_b(:,1:5) ; MESA_sys_myo_b(:,1:5)]; %from the PDM
trainingdata(:,7:11) = [DETERMINE_dia_myo_b(:,1:5) ; MESA_dia_myo_b(:,1:5)];
trainingdata(:,12:21) = [DETERMINE_dia_sys_myo_b(:,1:10) ; MESA_dia_sys_myo_b(:,1:10)];
%sphericity measure
trainingdata(:,22) = [DETERMINE_dia_endo_AVratios ; MESA_dia_endo_AVratios ];
trainingdata(:,23) = [DETERMINE_sys_endo_AVratios ; MESA_sys_endo_AVratios ];
% mean thickness
trainingdata(:,24) = [diaMeanT(data(1).DETERMINE_indices)' ; diaMeanT(data(1).MESA_indices)'];
trainingdata(:,25) = [sysMeanT(data(1).DETERMINE_indices)' ; sysMeanT(data(1).MESA_indices)'];
% mean(sys thickness) - mean(dia thickness) e.g. change in mean thickness from diastole to systole
trainingdata(:,26) = [deltaMeanT(data(1).DETERMINE_indices)' ; deltaMeanT(data(1).MESA_indices)'];
%most common thickness
trainingdata(:,27) = [diaModeT(data(1).DETERMINE_indices)' ; diaModeT(data(1).MESA_indices)'];
trainingdata(:,28) = [sysModeT(data(1).DETERMINE_indices)' ; sysModeT(data(1).MESA_indices)'];
% median thickness
trainingdata(:,29) = [sysMedianT(data(1).DETERMINE_indices)' ; sysMedianT(data(1).MESA_indices)'];
%thickness variance for diastole and systole seperately
trainingdata(:,30) = [DETERMINEdiaMyoThicknessVariance' ; MESAdiaMyoThicknessVariance' ];
trainingdata(:,31) = [DETERMINEsysMyoThicknessVariance' ; MESAsysMyoThicknessVariance' ];
% var(sys thickness - dia thickness)  e.g. variance in the difference between diastole and systole
trainingdata(:,32) = [DETERMINEmyoTchangesVar' ;MESAmyoTchangesVar'];
% % var([ dia thickness ; sys thickness])
% trainingdata(:,33) = [DETERMINE_DiaSys_dEPI2ENDO_vars' ; MESA_DiaSys_dEPI2ENDO_vars' ];
%% training data labels
trainingdataLabels = [
    'EF                ' ;...
    'sys myo b1        ' ;...
    'sys myo b2        ' ;...
    'sys myo b3        ' ;...
    'sys myo b4        ' ;...
    'sys myo b5        ' ;...
    'dia myo b1        ' ;...
    'dia myo b2        ' ;...
    'dia myo b3        ' ;...
    'dia myo b4        ' ;...
    'dia myo b5        ' ;...
    'dia sys myo b1    ' ;...
    'dia sys myo b2    ' ;...
    'dia sys myo b3    ' ;...
    'dia sys myo b4    ' ;...
    'dia sys myo b5    ' ;...
    'dia sys myo b6    ' ;...
    'dia sys myo b7    ' ;...
    'dia sys myo b8    ' ;...
    'dia sys myo b9    ' ;...
    'dia sys myo b10   ' ;...
    'dia endo AV ratio ' ;...
    'sys endo AV ratio ' ;...
    'diaMeanT          ' ;...
    'sysMeanT          ' ;...
    'deltameanT        ' ;...
    'diaModeT          ' ;...
    'sysModeT          ' ;...
    'sysMedianT        ' ;...
    'diaMyoThicknessVar' ;...
    'sysMyoThicknessVar' ;...
    'myoTchangesVar    ' ] ;

%% class names
names = char(200,1);
names(1:100,1) = 'd';
names(101:200,1) = 'm';

%% visualise training data
bestParam3 = 1 %EF
b = [bestParam1(1);bestParam2(1);bestParam3]';
figure
hold on
%histogram(trainingdata(1:100,b))
plot3D(trainingdata(1:100,b),'o')
plot3D(trainingdata(101:200,b),'o')
hold off
legend ' DETERMINE' 'MESA'
xlabel ([ num2str(bestParam1(1)) ])
ylabel ([ num2str(bestParam2(1)) ])
zlabel ([ num2str(bestParam3) ' (EF) '])

%% svmtrain (will be removed from later versions of matlab)
svmStruct = svmtrain(trainingdata(:,[bestParam1(1);bestParam2(1)]),names,'ShowPlot', true, 'kernel_function', 'linear');

%% SVM - 3 parameters
% %  which is the best trio of b values for classification?
tic
SVMclassification_rates = zeros((size(param1,2)),(size(param2,2)));
for param1 = 2:33
    for param2 = 29
        for   param3 = 1 % ejection fractions
            tic
            param1
            param2
            % SVMModel = fitcsvm(trainingdata, group, 'Standardize', true)
            SVMModel = fitcsvm(trainingdata(:,[param1;param2]), names, 'KernelFunction', 'linear' );
            %SVMModel = fitcsvm(trainingdata(:,:), names, 'KernelFunction', 'linear' );
            
            
            % cross-validation of the SVM
            %             CVSVMModel = crossval(SVMModel);
            %             misclassification_rate = kfoldLoss(CVSVMModel);
            %             %         classification_rates(param1,param2) = 1 - misclassification_rate
            %             SVMclassification_rates(param1, param2) = 1 - misclassification_rate
            %
        end
    end
end

toc
%% Linear discriminant analysis (LDA)

%% LDA - 2 parameters

LDAcross_validated_classification_rates =  zeros(size(trainingdata,2),size(trainingdata,2));
for param1 = 1:33
    param1
    for param2 = 1
        param2
        obj = fitcdiscr(trainingdata(:,[param1;param2]), names);
        % obj = fitcdiscr(trainingdata(:,:), names); %using all 12 gives 94% classified
        resuberror = resubLoss(obj); %proportion of misclassified observations
        LDAclassification_rates(param1,param2) = 1 - resuberror;
        
        cvmodel = crossval(obj,'leaveout','on');
        cverror = kfoldLoss(cvmodel);
        LDAcross_validated_classification_rates(param1, param2) = 1 - cverror;
        
        
    end
end

max_rate = max(max(LDAcross_validated_classification_rates))

figure
% subplot 121
plot([0 36],[max_rate max_rate])
legend (['max classification = ' num2str(max_rate)])
hold on
plot([0 36],[0.87 0.87])
bar( LDAcross_validated_classification_rates(:,1))
ylabel 'classification success'
ylim ([0.8 1])
xlim ([1 36])
set(gca,'XTick',[1:1:36])
% set(gca,'XTickLabel',trainingdataLabels)
xlabel 'parameter'
set(gca,'yMinorTick','on')
title 'EF + ...'
% subplot 122
% histogram(LDAclassification_rates)
% xlabel

[bestParam1,bestParam2] = find(LDAcross_validated_classification_rates==max(max(LDAcross_validated_classification_rates)))
[best1] = trainingdataLabels(bestParam1,:)
[best2] = trainingdataLabels(bestParam2,:)
%% LDA - 3 parameters
LDAcross_validated_classification_rates =  zeros(size(param1,2),size(param2,2));
for param1 = 1:33
    param1
    for param2 = 27
        for param3 = 1 % param3 is EF
            tic
            obj = fitcdiscr(trainingdata(:,[param1, param2 param3]), names);
            
            resuberror = resubLoss(obj); %proportion of misclassified observations
            LDAclassification_rates(param1, param2) = 1 - resuberror;
            cvmodel = crossval(obj,'leaveout','on');
            cverror = kfoldLoss(cvmodel);
            LDAcross_validated_classification_rates(param1, param2) = 1 - cverror;
            
        end
    end
end

hold on
% subplot 121
plot([0 size(trainingdata,2)],[max(max(LDAcross_validated_classification_rates)) max(max(LDAcross_validated_classification_rates))])
legend (['max classification = ' num2str(max(max(LDAcross_validated_classification_rates)))])
plot([0 size(trainingdata,2)],[0.87 0.87])
plot(LDAcross_validated_classification_rates(:,param2),'o')
ylabel 'classification success'
ylim ([0.8 1])
xlim ([1 size(trainingdata,2)])
xlabel 'training data column (parameter index)'
set(gca,'yMinorTick','on')
title 'EF + systolic endocardium area to volume ratio + ...'
% subplot 122

%parameters that give best classification when used with EF
max(max(LDAcross_validated_classification_rates))
[bestParam1,bestParam2] = find(LDAcross_validated_classification_rates==max(max(LDAcross_validated_classification_rates)))
[best1] = trainingdataLabels(bestParam1,:)
[best2] = trainingdataLabels(bestParam2,:)
% figure
% subplot 121
% bar3(LDAclassification_rates(:,:))
% zlim ([0.8 1])
% subplot 122
% imagesc((LDAclassification_rates(:,:)))
% title 'LDA - 3 parameters'
% xlabel 'training data column'
% ylabel 'training data column'
% colorbar

% LDAclassification_rates = table(LDAclassification_rates(:,:))


LDAcross_validated_classification_rates =  zeros(size(param1,2),size(param2,2));
for param4 = 1:33
    for param1 = 32
    param4
    for param2 = 27
        for param3 = 1 % param3 is EF
            tic
            obj = fitcdiscr(trainingdata(:,[param4, param1, param2 param3]), names);
            
            resuberror = resubLoss(obj); %proportion of misclassified observations
            LDAclassification_rates(param4) = 1 - resuberror;
            cvmodel = crossval(obj,'leaveout','on');
            cverror = kfoldLoss(cvmodel);
            LDAcross_validated_classification_rates(param4) = 1 - cverror;
            
        end
    end
end
end
hold on
% subplot 121
plot([0 size(trainingdata,2)],[max(max(LDAcross_validated_classification_rates)) max(max(LDAcross_validated_classification_rates))])
legend (['max classification = ' num2str(max(max(LDAcross_validated_classification_rates)))])
plot([0 size(trainingdata,2)],[0.87 0.87])
plot(LDAcross_validated_classification_rates(1:param4),'o')
ylabel 'classification success'
ylim ([0.8 1])
xlim ([1 size(trainingdata,2)])
xlabel 'training data column (parameter index)'
set(gca,'yMinorTick','on')
title 'EF + systolic endocardium area to volume ratio + ...'
% subplot 122
%% LDA - confusion matrix shows performance of classifier.
% row 1: DETERMINE, row 2: MESA
R = confusionmat(obj.Y,resubPredict(obj))

%% plot LDA example
sysMyob1 = trainingdata(:,1);
EF = trainingdata(:,11);
figure
h1 = gscatter(sysMyob1,EF,names,'krb','ov^',[],'off');
hold on
h1(1).LineWidth = 2;
h1(2).LineWidth = 2;
K = obj.Coeffs(1,2).Const;
L = obj.Coeffs(1,2).Linear;
f = @(sysmyob1,EF)K + L(1)*sysmyob1 + L(2)*EF;
h2 = ezplot(f,[1100 2000 -1300 300]);
h2.Color = 'g';
h2.LineWidth = 2;
legend 'DETERMINE' 'MESA' 'LDA classification boundary' 'Location' 'best'
hold off

%% visualise classification tables (SVM and LDA)
figure
subplot 121
% surf(classification_rates)
imagesc(SVMclassification_rates)
title ' classification rates'
% xlim ([1 ; 5]')
% ylim ([1 ; 5])
xlabel 'b'
ylabel 'b'
caxis([.5, .9])
colorbar

subplot 122
% surf(classification_rates_withEF(:,:,6))
imagesc(LDAclassification_rates_withEF(:,:,6))
title ' classification rates with EF'
xlabel 'b'
ylabel 'b'
caxis([.8, .9])
colorbar


%% kmeans
idx = kmeans(trainingdata(:,:),2);
% plot(names,idx)
% histogram(names);
% hold on
% histogram(idx);
clear SCORE
for i = 1:200
    if names(i) == idx(i)
        SCORE(i) = 1;
    end
end
classification_rate = (sum(SCORE)/200)

%% Decision tree
% rng(1);
% tree = fitctree(trainingdata(:,:), names,'CrossVal','on')
% numBranches = @(x)sum(x.IsBranch);
% mdlDefaultNumSplits = cellfun(numBranches, tree.Trained);

% figure;
% histogram(mdlDefaultNumSplits)

% view(tree.Trained{1},'Mode','graph')
%
% E = kfoldLoss(tree)
% classification_success = 1 - E


%random forest
nTrees =300
BaggedEnsemble = TreeBagger(nTrees, trainingdata(:,:),names,'method','classification','OOBPred','On') %'NVarToSample',100)

oobErrorBaggedEnsemble = oobError(BaggedEnsemble);
figure
plot(oobErrorBaggedEnsemble)
xlabel 'Number of grown trees';
ylabel 'Out-of-bag classification error';

% 1 - oobErrorBaggedEnsemble(nTrees,1)
%
% err = error(BaggedEnsemble,trainingdata(:,:),names,'mode','cumulative')
% plot(err)
% hold on
% plot([0 100] , [mean(err) mean(err)],'r')
% legend 'individual error' 'mean error'
% 1 - mean(err)
% hold off
% histogram(err)
% vals = crossval(BaggedEnsemble, trainingdata(:,:))
%% performance curves - comparing classifiers
% [X,Y] = perfcurve(names,scores,posclass)

% %Comparing classifiers
% resp(1:200,1:2) = strcmp(names,'d'); % resp = 1, if Y = 'b', or 0 if Y = 'g'
% pred = trainingdata(:,1:2);
% mdl = fitglm(pred,resp,'Distribution','binomial','Link','logit');
% score_log = mdl.Fitted.Probability; % Probability estimates
% [X,Y,T,AUC] = perfcurve(names,score_svm(:,mdlSVM.ClassNames),'true');
% AUC
% resp = strcmp(names(:,1:2),'b');
% pred = trainingdata(:,1:2);
% mdl = fitglm(pred,resp,'Distribution','binomial','Link','logit');
% score_log = mdl.Fitted.Probability;
% [Xsvm,Ysvm,Tsvm,AUCsvm] = perfcurve(resp,score_svm(:,mdlSVM.ClassNames),'true');

% pred = trainingdata(:,1:2);
%
% for i = 1:200
%     if names(i,1) == 'm'
%         resp(i,1) = 1
%     else resp(i,1) = 0
%     end
% end
%
% SVMModel_2 = fitPosterior(SVMModel_2);
% SVMModel_1 = fitPosterior(SVMModel_1);
% [~,scores2] = resubPredict(SVMModel_2);
% [~,scores1] = resubPredict(SVMModel_1);
% [x1,y1,~,auc1] = perfcurve(names,scores1(:,2),1);
% [x2,y2,~,auc2] = perfcurve(resp,scores2(:,2),1);