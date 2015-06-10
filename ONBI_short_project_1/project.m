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

%% histogram of SUM(systole thickness)/SUM(diastole thickness)
figure
histogram(MyoThicknessSysDiaRatio(data(1).MESA_indices))
hold on
histogram(MyoThicknessSysDiaRatio(data(1).DETERMINE_indices))
legend 'MESA' 'DETERMINE'
title 'sum(systole thickness)/sum(diastole thickness)'
xlabel 'sum(systole thickness)/sum(diastole thickness)'
ylabel 'frequency'

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
ylim  ([0 10e5]);
for i = [data(1).DETERMINE_indices ; data(1).MESA_indices]'
    jackVolumes1(i,1) = data(i).diastolic.myoVolume;
    jackVolumes2(i,1) = data(i).systolic.myoVolume;
end
%      jackVolumes1(~any( jackVolumes1,2),:) = [];
%      jackVolumes2(~any( jackVolumes2,2), : ) = [];
indices = [data(1).DETERMINE_indices ; data(1).MESA_indices]
hold on

x = [0 10e5]
plot(x, x)
plot(jackVolumes1([data(1).DETERMINE_indices ; data(1).MESA_indices])', ErnestoVolumes.data(:,4)', '.');
plot(jackVolumes2([data(1).DETERMINE_indices ; data(1).MESA_indices])', ErnestoVolumes.data(:,3)', '.');
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
DETERMINE_sys_epi_AVratios = cell2mat({data(:).DETERMINE_sys_epi_AVratio})' ;
DETERMINE_dia_epi_AVratios = cell2mat({data(:).DETERMINE_dia_epi_AVratio})' ;
MESA_sys_epi_AVratios = cell2mat({data(:).MESA_sys_epi_AVratio})' ;
MESA_dia_epi_AVratios = cell2mat({data(:).MESA_dia_epi_AVratio})' ;
figure
subplot 221
histogram(DETERMINE_dia_endo_AVratios)
hold on
histogram(MESA_dia_endo_AVratios)
title 'dia endo sphericity'
legend 'DETERMINE' 'MESA'
subplot 222
histogram(DETERMINE_dia_epi_AVratios)
hold on
histogram(MESA_dia_epi_AVratios)
title 'dia epi sphericity'
legend 'DETERMINE' 'MESA'
subplot 223
histogram(DETERMINE_sys_endo_AVratios)
hold on
histogram(MESA_sys_endo_AVratios)
title 'sys endo sphericity'
legend 'DETERMINE' 'MESA'
subplot 224
histogram(DETERMINE_sys_epi_AVratios)
hold on
histogram(MESA_sys_epi_AVratios)
title 'sys epi sphericity'
legend 'DETERMINE' 'MESA'

%% visualise mesh area stats
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
%% PCA - find eigenvectors of covariance matrix 
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
nEigenvalues = 15;
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
disp('start choosing principle eigenvectors')
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
disp('finished finding principle eigenvectors')
%% Find model parameter values (b) , using +/- sqrt(eigenvalue) as b
% b can be increased later.
disp('start finding max and min b values')
% systolic, endocardium

dia_myo_max_b = sqrt(principle_dia_myo_eigenvalues);

sys_myo_max_b = sqrt(principle_sys_myo_eigenvalues);

dia_sys_myo_max_b = sqrt(principle_dia_sys_myo_eigenvalues);
disp('finished finding max and min b values')
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
for b = 1:nEigenvalues
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
run('visualise_eigenmode_contributions.m')
%% visualise modes - moving
run('visualise_eigenmode_variations_moving.m')
%% visualise modes - static
run('visualise_eigenmode_variations_static.m')
%% play mode of variation movie
figure
title ([ 'mode',  num2str(visualisingMode) ])
axis off
for n = 1:4
    movie(F)
    
end


%% Bullseye plots
nStd = 3
% diaPts = [1:2178];
% sysPts = [2179:2*2178];
cmin = 0
cmax = 25
run('bullseye_plots.m')
visualisingMode1 
visualisingMode2 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% find myocardium thickness eigenmodes
run('find_myocardium_thickness_eigenmodes.m') ;%plots thicknesses and eigenmode contributions
visualisingMode1 
visualisingMode2 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CLASSIFICATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SVM - set training data and class names

run('set_training_data.m')

%% visualise b values 
b = [1;4;7];
run('visualise_b_values__histo_2D_3D.m')
%% visualise training data
run('visualise_training_data.m')

%% svmtrain (will be removed from later versions of matlab)
svmStruct = svmtrain(trainingdata(:,[bestParam1(1);bestParam2(1)]),names,'ShowPlot', true, 'kernel_function', 'linear');

%% SVM - 3 parameters
% %  which is the best trio of b values for classification?
run('svm_classification_3_params.m')

%% Linear discriminant analysis (LDA)
%% LDA - 2 parameters
run('lda_classification_2_params.m')

%% LDA - 3 parameters
run('lda_classification_3_params.m')

%% LDA - 4 parameters
run('lda_classification_4_params.m')

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

tic
% cvpart = cvpartition(names,'holdout',0.3);
% Xtrain = trainingdata(training(cvpart),:);
% Ytrain = names(training(cvpart),:);
% Xtest = trainingdata(test(cvpart),:);
% Ytest = names(test(cvpart),:);
paramIndices = [60:67];
trainingdataLabels(paramIndices,:)
rng(1);
BaggedEnsemble = TreeBagger(300, trainingdata(:,paramIndices),names,'method','classification','OOBPred','On','oobvarimp','on') %'NVarToSample',100)

% BaggedEnsemble = fitensemble(Xtrain,Ytrain,'Bag',200,'Tree','Type','Classification')
% oobErrorBaggedEnsemble = oobError(BaggedEnsemble);
%
% random forest - estimate feature importance
figure
% subplot 121
% "For any variable, the measure is the increase in prediction error if the
% values of that variable are permuted across the out-of-bag observations"
% "By default, NVarToSample is equal to the square root of the total number
% of variables for classification"
bar(BaggedEnsemble.OOBPermutedVarDeltaError);
sorted_importance = sort(BaggedEnsemble.OOBPermutedVarDeltaError,'descend');
find(BaggedEnsemble.OOBPermutedVarDeltaError == sorted_importance(1))
% important1 = trainingdataLabels(find(BaggedEnsemble.OOBPermutedVarDeltaError == sorted_importance(1)),:)
find(BaggedEnsemble.OOBPermutedVarDeltaError == sorted_importance(2))
% important2 = trainingdataLabels(find(BaggedEnsemble.OOBPermutedVarDeltaError == sorted_importance(2)),:)
find(BaggedEnsemble.OOBPermutedVarDeltaError == sorted_importance(3))
% important3 = trainingdataLabels(find(BaggedEnsemble.OOBPermutedVarDeltaError == sorted_importance(3)),:)
set(gca,'XTick',[1:size(BaggedEnsemble.OOBPermutedVarDeltaError,2)])
set(gca,'XTickLabel',trainingdataLabels(paramIndices,:))
xlabel('Feature index')
ylabel('out-of-bag feature importance')
% subplot 122
% hold on
% plot([0 300], [0.07 0.07])
% plot(oobErrorBaggedEnsemble)
% xlabel 'Number of grown trees';
% ylabel 'Out-of-bag classification error';
% legend ([num2str(0.07)])
toc
%%
tic
paramIndices = [1:68];
% for n = 1:2
%     paramIndices = paramIndices + (n-1)*33
% n

figure('name','random forest classification cross-validation');
% plot(loss(BaggedEnsemble,Xtest,Ytest,'mode','cumulative'));
xlabel('Number of trees');
ylabel('classification error');
hold on
%
%K-FOLD
%original sample is randomly partitioned into k equal sized subsamples
% k-1 of the subsamples are used as training data.
%cross-validation repeated k times
rng(1);
k = 5
kfold_cross_validatedBE = fitensemble(trainingdata(:,paramIndices),names,'Bag',300,'Tree','type','classification','kfold',k);
kfoldlossplot = kfoldLoss(kfold_cross_validatedBE,'mode','cumulative')
plot(kfoldlossplot,'g.'); %kfold =5

%  
%LEAVE-ONE-OUT
% when 'leaveout' is set to 'on' , kfold = number of observations 
rng(1);
leaveout_cross_validatedBE = fitensemble(trainingdata(:,paramIndices),names,'Bag',300,'Tree','type','classification','leaveout','on');
kfoldlossLOOplot = kfoldLoss(leaveout_cross_validatedBE, 'mode','cumulative')
plot(kfoldlossLOOplot,'k-'); %kfold =200
% imp = predictorImportance(leaveout_cross_validatedBE)
% figure; bar(imp)
% end
toc
%%
% plot(oobLoss(BaggedEnsemble,'mode','cumulative'),'g.');
line = 0.13;
plot([0 300] ,[line line])
% legend(['kfold, k = ',num2str(k),])
legend 'kfold, k = 5' ' leave one out' 'error = 0.13'
% 'leave-one-out ', num2str(line)']); 

%%
% cvmodel = crossval(obj,'leaveout','on');
% cverror = kfoldLoss(cvmodel);
% LDAcross_validated_classification_rates(param1, param2) = 1 - cverror;
%         

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