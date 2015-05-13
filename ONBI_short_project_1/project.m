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
addpath C:\Users\jesu2687\Documents\MATLAB\ONBI_short_project_1\output-stacom-newcase\output-stacom-newcase
load('C:\Users\jesu2687\Documents\MATLAB\ONBI_short_project_1\project1_data.mat');
load('C:\Users\jesu2687\Documents\MATLAB\ONBI_short_project_1\labels.mat');

%store indices that allow us to indentify which study (DETERMINE or MESA)
%each case belongs to.
data(1).DETERMINE_indices = Labels(2:101,3);
data(1).MESA_indices = Labels(102:201,3);
data(1).MESA_indices(2) = 401 ; % replace SMM0001 with the replacement provided by the organisers (SMM0401
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  CALCULATE EJECTION FRACTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('started calculating ejection fractions')

[data] = calcEjectionFraction(data);

disp('finished calculating ejection fractions')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  CALCULATE ACCURACIES 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[data, accuracy, sensitivity, specificity] = calcAccuracy(data);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATE MYOCARDIUM VOLUMES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('calculating myocardium volumes')

%load struct containing the triangle file for the myo 'donut' shaped lid
load('C:\Users\jesu2687\Documents\MATLAB\ONBI_short_project_1\myoB_tri.mat') 

%calculate myocardium volumes
[data] = calcMyoVolume(data, myoB);

% store volumes for plotting
[data] = storeMyoVolumes(data);

disp('finished calculating myocardium volumes')

% plot systolic myocardium volumes
plotMyoVolumes(data)

% calculate "myocardium ejection fraction" myoEF
% myoSV = diastolic myovolume - systolic myovolume
% myoEF = (myoSV/(diastolic volume))*100
[data] = calcMyoEjectionFractions(data);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATE MESH AREAS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate triangle side lengths
%!!!!! NEED TO CORRECT THIS TO INCLUDE THE LIDS
for i = 1:400 
    for d = find(MESA_indices==i)
       
       MESA_dia_endo_sides = calcTriSides(data(i).diastolic.endo.tri, data(i).diastolic.endo.xyz); 
       MESA_dia_endo_areas(d,1) = calcTriMeshArea(MESA_dia_endo_sides);
          
       MESA_dia_epi_sides = calcTriSides(data(i).diastolic.epi.tri, data(i).diastolic.epi.xyz);
       MESA_dia_epi_areas(d,1) = calcTriMeshArea(MESA_dia_epi_sides);
       
       MESA_sys_endo_sides = calcTriSides(data(i).systolic.endo.tri, data(i).systolic.endo.xyz);
       MESA_sys_endo_areas(d,1) = calcTriMeshArea(MESA_sys_endo_sides);
       
       MESA_sys_epi_sides = calcTriSides(data(i).systolic.epi.tri, data(i).systolic.epi.xyz);
       MESA_sys_epi_areas(d,1) = calcTriMeshArea(MESA_sys_epi_sides);
       
       MESA_dia_myo_sides = calcTriSides(data(i).diastolic.myo.tri, data(i).diastolic.myo.xyz);
       MESA_dia_myo_areas(d,1) = calcTriMeshArea(MESA_dia_myo_sides);
       
       MESA_sys_myo_sides = calcTriSides(data(i).systolic.myo.tri, data(i).systolic.myo.xyz);
       MESA_sys_myo_areas(d,1) = calcTriMeshArea(MESA_sys_myo_sides);
    end
end
for i = 1:400
    for d = find(DETERMINE_indices==i)
       DETERMINE_dia_endo_sides = calcTriSides(data(i).diastolic.endo.tri, data(i).diastolic.endo.xyz); 
       DETERMINE_dia_endo_areas(d,1) = calcTriMeshArea(DETERMINE_dia_endo_sides);
          
       DETERMINE_dia_epi_sides = calcTriSides(data(i).diastolic.epi.tri, data(i).diastolic.epi.xyz);
       DETERMINE_dia_epi_areas(d,1) = calcTriMeshArea(DETERMINE_dia_epi_sides);
       
       DETERMINE_sys_endo_sides = calcTriSides(data(i).systolic.endo.tri, data(i).systolic.endo.xyz);
       DETERMINE_sys_endo_areas(d,1) = calcTriMeshArea(DETERMINE_sys_endo_sides);
       
       DETERMINE_sys_epi_sides = calcTriSides(data(i).systolic.epi.tri, data(i).systolic.epi.xyz);
       DETERMINE_sys_epi_areas(d,1) = calcTriMeshArea(DETERMINE_sys_epi_sides);
       
       DETERMINE_dia_myo_sides = calcTriSides(data(i).diastolic.myo.tri, data(i).diastolic.myo.xyz);
       DETERMINE_dia_myo_areas(d,1) = calcTriMeshArea(DETERMINE_dia_myo_sides);
       
       DETERMINE_sys_myo_sides = calcTriSides(data(i).systolic.myo.tri, data(i).systolic.myo.xyz);
       DETERMINE_sys_myo_areas(d,1) = calcTriMeshArea(DETERMINE_sys_myo_sides);
    end
end

%% plot surface areas
figure;
hold on
nbins = 25;
histogram(DETERMINE_dia_endo_areas,nbins)
histogram(MESA_dia_endo_areas,nbins)
legend 'DETERMINE' ' MESA'
title 'diastolic endocardium surface areas'
xlabel 'surface areas (mm^2)'
ylabel 'frequency'
print('compare diastolic endocardium surface areas','-dpng')

figure;
hold on
nbins = 25;
histogram(DETERMINE_dia_epi_areas,nbins)
histogram(MESA_dia_epi_areas,nbins)
legend 'DETERMINE' ' MESA'
title 'diastolic epicardium surface areas'
xlabel 'surface areas (mm^2)'
ylabel 'frequency'
print('compare diastolic epicardium surface areas','-dpng')

figure;
hold on
nbins = 25;
histogram(DETERMINE_sys_endo_areas,nbins)
histogram(MESA_sys_endo_areas,nbins)
legend 'DETERMINE' ' MESA'
title 'systolic endocardium surface areas'
xlabel 'surface areas (mm^2)'
ylabel 'frequency'
print('compare systolic endocardium surface areas','-dpng')


figure;
hold on
nbins = 25;
histogram(DETERMINE_sys_epi_areas,nbins)
histogram(MESA_sys_epi_areas,nbins)
legend 'DETERMINE' ' MESA'
title 'systolic epicardium surface areas'
xlabel 'surface areas (mm^2)'
ylabel 'frequency'
print('compare systolic epicardium surface areas','-dpng')

figure;
hold on
nbins = 25;
histogram(DETERMINE_dia_myo_areas,nbins)
histogram(MESA_dia_myo_areas,nbins)
legend 'DETERMINE' ' MESA'
title 'diastolic myocardium surface areas'
xlabel 'surface areas (mm^2)'
ylabel 'frequency'
print('compare diastolic myocardium surface areas','-dpng')

figure;
hold on
nbins = 25;
histogram(DETERMINE_sys_myo_areas,nbins)
histogram(MESA_sys_myo_areas,nbins)
legend 'DETERMINE' ' MESA'
title 'systolic myocardium surface areas'
xlabel 'surface areas (mm^2)'
ylabel 'frequency'
print('compare systolic myocardium surface areas','-dpng')

%% calculate surface area to volume ratios
for i = 1:100 %100 cases in each class
    DETERMINE_dia_endo_AVratio(i,1) = DETERMINE_dia_endo_areas(i,1)/DETERMINE_diastolic_endoVolumes(i,1);
    DETERMINE_dia_epi_AVratio(i,1) = DETERMINE_dia_epi_areas(i,1)/DETERMINE_diastolic_epiVolumes(i,1);
    DETERMINE_sys_endo_AVratio(i,1) = DETERMINE_sys_endo_areas(i,1)/DETERMINE_systolic_endoVolumes(i,1);
    DETERMINE_sys_epi_AVratio(i,1) = DETERMINE_sys_epi_areas(i,1)/DETERMINE_systolic_epiVolumes(i,1);
    
    MESA_dia_endo_AVratio(i,1) = MESA_dia_endo_areas(i,1)/MESA_diastolic_endoVolumes(i,1);
    MESA_dia_epi_AVratio(i,1) = MESA_dia_epi_areas(i,1)/MESA_diastolic_epiVolumes(i,1);
    MESA_sys_endo_AVratio(i,1) = MESA_sys_endo_areas(i,1)/MESA_systolic_endoVolumes(i,1);
    MESA_sys_epi_AVratio(i,1) = MESA_sys_epi_areas(i,1)/MESA_systolic_epiVolumes(i,1);
    
    MESA_dia_myo_AVratio(i,1) = MESA_dia_myo_areas(i,1)/MESA_diastolic_myovolumes(i,1);
    MESA_sys_myo_AVratio(i,1) = MESA_sys_myo_areas(i,1)/MESA_systolic_myovolumes(i,1);
    
    DETERMINE_dia_myo_AVratio(i,1) = DETERMINE_dia_myo_areas(i,1)/DETERMINE_diastolic_myovolumes(i,1);
    DETERMINE_sys_myo_AVratio(i,1) = DETERMINE_sys_myo_areas(i,1)/DETERMINE_systolic_myovolumes(i,1);
    
    
end

figure;
hold on
nbins = 25;
histogram(DETERMINE_dia_endo_AVratio,nbins)
histogram(MESA_dia_endo_AVratio,nbins)
legend 'DETERMINE' ' MESA'
title 'diastolic endocardium surface area to volume ratio'
xlabel 'area/volume (mm^-1)'
ylabel 'frequency'
print('compare diastolic endocardium surface areas to volume ratios','-dpng')

figure;
hold on
nbins = 25;
histogram(DETERMINE_sys_epi_AVratio,nbins)
histogram(MESA_sys_epi_AVratio,nbins)
legend 'DETERMINE' ' MESA'
title 'systolic endocardium surface area to volume ratio'
xlabel 'area/volume (mm^-1)'
ylabel 'frequency'
print('compare systolic endocardium surface area to volume ratios','-dpng')

figure;
hold on
nbins = 25;
histogram(DETERMINE_sys_endo_AVratio,nbins)
histogram(MESA_sys_endo_AVratio,nbins)
legend 'DETERMINE' ' MESA'
title 'systolic endocardium surface area to volume ratio'
xlabel 'area/volume (mm^-1)'
ylabel 'frequency'
print('compare systolic endocardium surface area to volume ratios','-dpng')

figure;
hold on
nbins = 25;
histogram(DETERMINE_sys_epi_AVratio,nbins)
histogram(MESA_sys_epi_AVratio,nbins)
legend 'DETERMINE' ' MESA'
title 'systolic epicardium surface area to volume ratio'
xlabel 'area/volume (mm^-1)'
ylabel 'frequency'
print('compare systolic epicardium surface area to volume ratios','-dpng')

figure;
hold on
nbins = 25;
histogram(DETERMINE_dia_myo_AVratio,nbins)
histogram(MESA_dia_myo_AVratio,nbins)
legend 'DETERMINE' ' MESA'
title 'diastolic myocardium surface area to volume ratio'
xlabel 'area/volume (mm^-1)'
ylabel 'frequency'
print('compare diastolic myocardium surface areas to volume ratios','-dpng')

%% calculate total areas(heron's forumla)
%!!!!! NEED TO CORRECT THIS TO INCLUDE THE LIDS
% endo_area = calcTriMeshArea(endo_sides);
% epi_area = calcTriMeshArea(epi_sides);
% %% calculate coordinates of the centroids
% %!!!!! NEED TO CORRECT THIS...ALREADY DONE BY 'calcVolume'?
% endo_centroid = mean(endo);
% epi_centroid = mean(epi);
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SUPPORT VECTOR MACHINE (SVM) CLASSIFICATION - VOLUMES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot data to perform svm on
hold on
% plot3(DETERMINE_EF, DETERMINE_systolic_endoVolumes,DETERMINE_diastolic_endoVolumes,'o')
% plot3(MESA_EF, MESA_systolic_endoVolumes,MESA_diastolic_endoVolumes,'o')
plot(DETERMINE_EF, DETERMINE_systolic_endoVolumes,'o')
plot(MESA_EF, MESA_systolic_endoVolumes,'o')
hold off

%% train SVM
trainingdata(:,1) = [DETERMINE_EF ; MESA_EF];
trainingdata(:,2) = [DETERMINE_systolic_endoVolumes ; MESA_systolic_endoVolumes];
%must be column vector, with each row corresponding to the value of the 
%corresponding row in trainingdata
classNames(1:100,1) = 'd';
classNames(101:200,1) = 'm';
% names(1:100,2) = 1;
% names(101:200,2) = 2;

% svmStruct = svmtrain(trainingdata,group,'ShowPlot', true);

% SVMModel = fitcsvm(trainingdata, group, 'Standardize', true)
SVMModel = fitcsvm(trainingdata, classNames,'KernelScale','auto','Standardize',true); %'ClassNames',{'d','m'});
classOrder = SVMModel.ClassNames
figure
gscatter(trainingdata(:,1),trainingdata(:,2),classNames)
hold on
sv = SVMModel.SupportVectors;
plot(sv(:,1),sv(:,2),'ko','MarkerSize',10)
legend('DETERMINE','MESA','Support Vector')
title 'SVM for EF and sys endo volume'
hold off

%% cross-validate the SVM classifier
CVSVMModel = crossval(SVMModel_1);
%calculate classification error
misclass = kfoldLoss(CVSVMModel_1);
misclassification_rate = misclass

%% SVM prediction
[label,score] = predict(SVMModel,newX);

%% SVM: EF and myo EF
% hold on
% % plot3(DETERMINE_EF, DETERMINE_systolic_endoVolumes,DETERMINE_diastolic_endoVolumes,'o')
% % plot3(MESA_EF, MESA_systolic_endoVolumes,MESA_diastolic_endoVolumes,'o')
% plot(DETERMINE_EF, DETERMINE_myoEF,'o')
% plot(MESA_EF, MESA_myoEF,'o')

%train the SVM
trainingdata(:,1) = [DETERMINE_EF ; MESA_EF];
trainingdata(:,2) = [DETERMINE_myoEF ; MESA_myoEF];
names = char(200,1);
names(1:100,1) = 'd';
names(101:200,1) = 'm';
% names(1:100,2) = 1;
% names(101:200,2) = 2;

% svmtrain will be removed from later versions of matlab.
% svmStruct = svmtrain(trainingdata,names,'ShowPlot', true);

% SVMModel = fitcsvm(trainingdata, group, 'Standardize', true)
SVMModel = fitcsvm(trainingdata, names, 'KernelFunction', 'linear' )
classOrder = SVMModel.ClassNames
figure
title 'SVM for EF and myoEF'
gscatter(trainingdata(:,1),trainingdata(:,2),names)
hold on
sv = SVMModel.SupportVectors;
plot(sv(:,1),sv(:,2),'ko','MarkerSize',10)
legend('DETERMINE','MESA','Support Vector')
% x = 10:80;
% plot(10:80, (0.005*x)+ -7.0407);
hold off


%% cross-validate the SVM
%"Determine the out-of-sample misclassification rate"
CVSVMModel = crossval(SVMModel);
misclassification_rate = kfoldLoss(CVSVMModel)


%% performance curves
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GENERAL PROCRUSTES ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is done by concatenating the endo and epi shape vectors to produce
% longer vectors (myo.xyz). Procrustes analysis is then performed on the
% two myo.xyz vectors. Finally, the transformed myo.xyz vectors are split
% to extract transformed endo and epi shape vectors.
disp('starting procrustes analysis')

%store dimensions of the shapes
[shape_nRows , shape_nCols] = size(data(2).diastolic.endo.xyz);

%use the initial data clouds as references (individual)
initial_reference_id = 2
dia_endo_reference = data(initial_reference_id).diastolic.endo.xyz;
dia_epi_reference = data(initial_reference_id).diastolic.epi.xyz;
sys_endo_reference = data(initial_reference_id).systolic.endo.xyz;
sys_epi_reference = data(initial_reference_id).systolic.epi.xyz;
%make a reference from the endo and epi surfaces (concatenate)
dia_myo_reference = [dia_endo_reference(:) ; dia_epi_reference(:)];
sys_myo_reference = [sys_endo_reference(:)  ; sys_epi_reference(:)];

dia_myo_reference = reshape(dia_myo_reference, [2*1089 3]);
sys_myo_reference = reshape(sys_myo_reference, [2*1089 3]);
%initialise the sums of the shapes as zero
% dia_endo_sum = zeros(size(dia_endo_reference));
% dia_epi_sum = zeros(size(dia_epi_reference));
% sys_endo_sum = zeros(size(sys_endo_reference));
% sys_epi_sum = zeros(size(sys_epi_reference));
% (whole shape)
dia_myo_sum = zeros(size(dia_myo_reference));
sys_myo_sum = zeros(size(sys_myo_reference));


% Think of endo and epi as one shape (concatenate).
for i = 1:401 %SMM001 has already been replaced by SMM401 in the script
data(i).diastolic.myo.xyz = [data(i).diastolic.endo.xyz ; data(i).diastolic.epi.xyz];
data(i).systolic.myo.xyz = [data(i).systolic.endo.xyz ; data(i).systolic.epi.xyz];
end

% MESA_indicesEDIT = [MESA_indices(1);MESA_indices(3:end)];
concatIndices = [data(1).MESA_indices ; data(1).DETERMINE_indices];
%p = number of times procrustes is performed.
for p = 1:3
procrustes_iteration = p
% iterate so that procrustes is performed on each case (400 patients).
for i = 1:400
%     i
%     for d = 1:find(concatIndices'==i)
% % transform endo and epi (individually)
% [data(i).diastolic.endo.procrustes_d, data(i).diastolic.endo.xyz] = procrustes(dia_endo_reference, data(i).diastolic.endo.xyz(:));
% [data(i).diastolic.epi.procrustes_d, data(i).diastolic.epi.xyz] = procrustes(dia_epi_reference, data(i).diastolic.epi.xyz(:));
% [data(i).systolic.endo.procrustes_d, data(i).systolic.endo.xyz] = procrustes(sys_endo_reference, data(i).systolic.endo.xyz(:));
% [data(i).systolic.epi.procrustes_d, data(i).systolic.epi.xyz] = procrustes(sys_epi_reference, data(i).systolic.epi.xyz(:));

% transform endo and epi as one shape vector
% [data(i).diastolic.myo.procrustes_d, data(i).diastolic.myo.xyz, dia_myo_transform] = procrustes(dia_myo_reference, data(i).diastolic.myo.xyz,'scaling',false,'reflection',false);
% [data(i).systolic.myo.procrustes_d, data(i).systolic.myo.xyz, sys_myo_transform] = procrustes(sys_myo_reference, data(i).systolic.myo.xyz,'scaling',false,'reflection',false);
[data(i).diastolic.myo.procrustes_d, diastolic_myo_shapes(i).xyz, dia_myo_transform(i)] = procrustes(dia_myo_reference, data(i).diastolic.myo.xyz,'scaling',false,'reflection',false);
[data(i).systolic.myo.procrustes_d, systolic_myo_shapes(i).xyz, sys_myo_transform(i)] = procrustes(sys_myo_reference, data(i).systolic.myo.xyz,'scaling',false,'reflection',false); %,'scaling',false); % with scaling SMM0362 goes much smaller than the others...


all_diastolic_myo_shapes(i,:) = diastolic_myo_shapes(i).xyz(:);
all_systolic_myo_shapes(i,:) = systolic_myo_shapes(i).xyz(:);

%sums for finding new means later on
% dia_endo_sum = data(i).diastolic.endo.xyz + dia_endo_sum;
% dia_epi_sum = data(i).diastolic.epi.xyz + dia_epi_sum;
% sys_endo_sum = data(i).systolic.endo.xyz + sys_endo_sum;
% sys_epi_sum = data(i).systolic.epi.xyz + sys_epi_sum;

% dia_myo_sum = data(i).diastolic.myo.xyz + dia_myo_sum;
% sys_myo_sum = data(i).systolic.myo.xyz + sys_myo_sum;
%      end
end

%calculate means
% dia_endo_mean = dia_endo_sum/400;
% dia_epi_mean = dia_epi_sum/400;
% sys_endo_mean = sys_endo_sum/400;
% sys_epi_mean = sys_epi_sum/400;
dia_myo_mean = diastolic_myo_shapes(1).xyz;
sys_myo_mean = systolic_myo_shapes(1).xyz;

for iS = 2:size(diastolic_myo_shapes,2)
    dia_myo_mean = dia_myo_mean + diastolic_myo_shapes(iS).xyz;
    sys_myo_mean = sys_myo_mean + systolic_myo_shapes(iS).xyz;
end

dia_myo_mean = dia_myo_mean./(size(diastolic_myo_shapes,2)-1);

sys_myo_mean = sys_myo_mean./(size(systolic_myo_shapes,2)-1);

%set means as the references for next round of procrustes
% dia_endo_reference = dia_endo_mean;
% dia_epi_reference = dia_epi_mean;
% sys_endo_reference = sys_endo_mean;
% sys_epi_reference = sys_epi_mean;
dia_myo_reference(i) = dia_myo_mean(i); 
sys_myo_reference(i) = sys_myo_mean(i);


end


% % reshape from vector to matrix.
% % individual endo and epi shapes.
% dia_endo_mean = reshape(dia_endo_mean,size(data(1).diastolic.endo.xyz));
% dia_epi_mean = reshape(dia_epi_mean,size(data(1).diastolic.endo.xyz));
% sys_endo_mean = reshape(sys_endo_mean,size(data(1).diastolic.endo.xyz));
% sys_epi_mean = reshape(sys_epi_mean,size(data(1).diastolic.endo.xyz));
dia_myo_reference = reshape(dia_myo_reference, [2*shape_nRows, shape_nCols]);
sys_myo_reference = reshape(sys_myo_reference, [2*shape_nRows, shape_nCols]);

% % concatenated endo and epi shapes.
% % split the long vectors into vectors that represent the endo and epi
% % shapes then reshape the resulting vectors to matrix form.
for i = 1:400
    

% d = find(concatIndices'==i);

% i
% reshape myo shapes from vector to matrix
% data(i).diastolic.myo.xyz = reshape(data(i).diastolic.myo.xyz,[2*shape_nRows, shape_nCols]);
% data(i).systolic.myo.xyz = reshape(data(i).systolic.myo.xyz,[2*shape_nRows, shape_nCols]);

% diastolic_myo_reshaped(i).xyz = reshape(diastolic_myo_shapes(:,i),size(dia_myo_reference));
% systolic_myo_reshaped(i).xyz = reshape(systolic_myo_shapes(:,i),size(sys_myo_reference));


% extract endo and epi shapes from full shape matrices (data(i).diastolic.myo.xyz)
data(i).diastolic.endo.xyz = data(i).diastolic.myo.xyz(1:shape_nRows, :);
data(i).diastolic.epi.xyz = data(i).diastolic.myo.xyz(shape_nRows+1:2*shape_nRows, :);
data(i).systolic.endo.xyz = data(i).systolic.myo.xyz(1:shape_nRows, :);
data(i).systolic.epi.xyz = data(i).systolic.myo.xyz(shape_nRows+1:2*shape_nRows, :);

end
disp('finished procrustes analysis')

%% Visualise procrustes

disp('visualise procrustes output')

plot3D(dia_myo_reference)
plot3D(sys_myo_reference)

% plotMESH = @(M,varargin)patch('vertices',M.xyz,'faces',M.tri,'edgecolor','k','facecolor','b',varargin{:});

%reference endocardium
%diastolic
figure
plot3(dia_endo_reference(:,1),dia_endo_reference(:,2),dia_endo_reference(:,3),'+')
axis equal; axis tight;
title 'reference endocardium (diastolic)'
%systolic
figure
plot3(sys_endo_reference(:,1),sys_endo_reference(:,2),sys_endo_reference(:,3),'+')
axis equal; axis tight;
title 'reference endocardium (diastolic)'

%reference epicardium
%distolic
figure
hold on
plot3(dia_epi_reference(:,1),dia_epi_reference(:,2),dia_epi_reference(:,3),'+')
axis equal; axis tight;
title 'reference epicardium (diastolic)'
%systolic
figure
plot3(sys_epi_reference(:,1),sys_epi_reference(:,2),sys_epi_reference(:,3),'+')
axis equal; axis tight;
title 'reference epicardium (systolic)'

%reference (whole LV shape from patient 1)
figure
plot3(dia_myo_reference(:,1), dia_myo_reference(:,2), dia_myo_reference(:,3),'o');
title 'reference LV'

%patient i endocardium and epicardium
%diastolic
figure
hold on
plot3(data(i).diastolic.endo.xyz(:,1), data(i).diastolic.endo.xyz(:,2), data(i).diastolic.endo.xyz(:,3),'o'); 
plot3(data(i).diastolic.epi.xyz(:,1), data(i).diastolic.epi.xyz(:,2), data(i).diastolic.epi.xyz(:,3),'o'); 
title 'patient i endo and epi, diastolic'
%systolic
figure
hold on
plot3(data(i).systolic.endo.xyz(:,1), data(i).systolic.endo.xyz(:,2), data(i).systolic.endo.xyz(:,3),'o'); 
plot3(data(i).systolic.epi.xyz(:,1), data(i).systolic.epi.xyz(:,2), data(i).systolic.epi.xyz(:,3),'o');
title 'patient i endo and epi, systolic'

%patient i whole LV shape
%diastolic
figure
hold on
plot3(data(i).diastolic.myo.xyz(:,1), data(i).diastolic.myo.xyz(:,2), data(i).diastolic.myo.xyz(:,3),'o'); 
title 'patient i, whole LV, diastolic'
%systolic

figure
for i = sort(MESA_indices')
    i
subplot 121
hold on
% plot3(dia_myo_reference(:,1), dia_myo_reference(:,2), dia_myo_reference(:,3),'go');
 plot3(data(i).systolic.endo.xyz(:,1), data(i).systolic.endo.xyz(:,2), data(i).systolic.endo.xyz(:,3),'.'); 
% plot3(data(i).systolic.epi.xyz(:,1), data(i).systolic.epi.xyz(:,2), data(i).systolic.epi.xyz(:,3),'.'); 
% plot3(data(i).systolic.myo.xyz(:,1), data(i).systolic.myo.xyz(:,2), data(i).systolic.myo.xyz(:,3),'.'); 

% plot3D(systolic_myo_reshaped(i).xyz)

title 'patient i, whole LV, systolic'
subplot 122
hold on
% plot3(dia_myo_reference(:,1), dia_myo_reference(:,2), dia_myo_reference(:,3),'go');
 plot3(data(i).diastolic.endo.xyz(:,1), data(i).diastolic.endo.xyz(:,2), data(i).diastolic.endo.xyz(:,3),'.'); 
% plot3(data(i).diastolic.epi.xyz(:,1), data(i).diastolic.epi.xyz(:,2), data(i).diastolic.epi.xyz(:,3),'.');
% plot3(data(i).diastolic.myo.xyz(:,1), data(i).diastolic.myo.xyz(:,2), data(i).diastolic.myo.xyz(:,3),'.'); 

% plot3D(diastolic_myo_reshaped(i).xyz)

title 'patient i, whole LV, systolic'
pause;

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SHAPE MODELING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('start shape modeling')
%% PCA
% see Cootes tutorial: http://personalpages.manchester.ac.uk/staff/timothy.f.cootes/Models/app_models.pdf

% %built in matlab method for PCA
% coeff = pca(data(1).diastolic.endo.xyz)
% figure
% hold on
% plot3(data(1).diastolic.endo.xyz(:,1), data(1).diastolic.endo.xyz(:,2), data(1).diastolic.endo.xyz(:,3),'o'); 
% plot3(coeff(:,1),coeff(:,2),coeff(:,3))

%% find mean shape and then use it to find the covariance matrix (of just MESA)
disp('started calculating covariance and mean shape')

dia_myo_cov_mat = cov(all_diastolic_myo_shapes);
sys_myo_cov_mat = cov(all_systolic_myo_shapes);
[dia_myo_principle_eigenvectors,dia_myo_principle_eigenvalues]=eigs(dia_myo_cov_mat,4);
[sys_myo_principle_eigenvectors,sys_myo_principle_eigenvalues]=eigs(sys_myo_cov_mat,4);
 
disp('finished calculating covariance and mean shape')
%% find eigenvectors of covariance matrix (of just MESA)
disp('started calculating eigenvectors of covariance matrix')

[dia_myo_principle_eigenvectors,dia_myo_principle_eigenvalues]=eigs(dia_myo_cov_mat,4);
[sys_myo_principle_eigenvectors,sys_myo_principle_eigenvalues]=eigs(sys_myo_cov_mat,4);

% % Each eigenvalue gives the variance of the data about the mean in the
% % direction of the corresponding eigenvector. Compute the total variance
%  totalVariance = sum(dia_endo_eigenvalues(:));
% proportion = 1;
% % %Choose the first t largest eigenvalues such that
%  p = proportion*totalVariance;
% a = dia_endo_eigenvalues/sum(dia_endo_eigenvalues(:));

disp('finished calculating eigenvectors of covariance matrix')
%% Find b values - using +/- 3*sqrt(eigenvalue) as b
% systolic, endocardium
dia_myo_min_b = - 3*sqrt(dia_myo_principle_eigenvalues);
dia_myo_max_b = 3*sqrt(dia_myo_principle_eigenvalues);
sys_myo_min_b = - 3*sqrt(sys_myo_principle_eigenvalues);
sys_myo_max_b = 3*sqrt(sys_myo_principle_eigenvalues);

%% (finding b by analysing each of the training shapes)
% for i = 1:400
%     dia_endo_b(i) = (dia_endo_principle_eigenvectors')*(data(i).diastolic.endo.xyz(:) - dia_endo_mean_shape);
%     sys_endo_b(i) = (sys_endo_principle_eigenvectors')*(data(i).systolic.endo.xyz(:) - sys_endo_mean_shape);
% end
% 
% % sys_endo_max_b = max(sys_endo_b);
% % sys_endo_min_b = min(sys_endo_b);
% 
% % plot(b)
% figure
% hold on
% histogram(dia_endo_b,20);
% histogram(sys_endo_b,20);
% legend 'diastolic b' 'systolic b'
% 
% dia_endo_max_b = std(dia_endo_b);
% dia_endo_min_b = -std(dia_endo_b);
% sys_endo_max_b = std(sys_endo_b);
% sys_endo_min_b = -std(sys_endo_b);

%% Calculate new shapes for visualisation
% mean shape + the sum of each eigenvector times it's value of b
%set of the principle eigenvectors to use (eigenvalue_range)
%(can equal the number_of_eigenvectors already extracted, or we can vary it for plotting)
nEigenmodes = 2;
for n = 1 : nEigenmodes
% dia_endo_new_shape_min(:,n) = dia_endo_mean_shape + dia_endo_principle_eigenvectors(:,eigenvalues(1,n)).*dia_endo_min_b(1,eigenvalues(1,n));
% dia_endo_new_shape_max(:,n) = dia_endo_mean_shape + dia_endo_principle_eigenvectors(:,eigenvalues(1,n)).*dia_endo_max_b(1,eigenvalues(1,n));
% sys_endo_new_shape_min(:,n) = sys_endo_mean_shape + sys_endo_principle_eigenvectors(:,eigenvalues(1,n)).*sys_endo_min_b(1,eigenvalues(1,n));
% sys_endo_new_shape_max(:,n) = sys_endo_mean_shape + sys_endo_principle_eigenvectors(:,eigenvalues(1,n)).*sys_endo_max_b(1,eigenvalues(1,n));
sys_myo_new_shape_min(:,n) = sys_myo_mean(:) + sys_myo_principle_eigenvectors(:,n)*sys_myo_min_b(n,n);
sys_myo_new_shape_max(:,n) = sys_myo_mean(:) + sys_myo_principle_eigenvectors(:,n)*sys_myo_max_b(n,n);
end
%% visualise eigenmodes of shape variation - points
% dia_endo_new_shape_min = reshape(dia_endo_new_shape_min(:,eigenmode), [1089 3]);
% dia_endo_new_shape_max = reshape(dia_endo_new_shape_max(:,eigenmode), [1089 3]);
% dia_endo_mean_shape = reshape(dia_endo_mean_shape(:,eigenmode), [1089 3]);
% sys_endo_new_shape_min = reshape(sys_endo_new_shape_min(:,eigenmode), [1089 3]);
% sys_endo_new_shape_max = reshape(sys_endo_new_shape_max(:,eigenmode), [1089 3]);
% sys_endo_mean_shape = reshape(sys_endo_mean_shape(:,eigenmode), [1089 3]);

sys_myo_new_shape_1_min = reshape(sys_myo_new_shape_min(:,1), [2178 3]);
sys_myo_new_shape_1_max = reshape(sys_myo_new_shape_max(:,1), [2178 3]);
sys_myo_new_shape_2_min = reshape(sys_myo_new_shape_min(:,2), [2178 3]);
sys_myo_new_shape_2_max = reshape(sys_myo_new_shape_max(:,2), [2178 3]);

figure
hold on
title 'PCA - eigenmode visualisation - diastolic endocardium'
plot3(dia_endo_mean_shape(:,1),dia_endo_mean_shape(:,2),dia_endo_mean_shape(:,3),'g.')
plot3(dia_endo_new_shape_max(:,1),dia_endo_new_shape_max(:,2),dia_endo_new_shape_max(:,3),'ro')
plot3(dia_endo_new_shape_min(:,1),dia_endo_new_shape_min(:,2),dia_endo_new_shape_min(:,3),'bo')
legend 'mean' 'max b' 'min b'

figure
title 'PCA - eigenmode visualisation - systolic endocardium'
hold on
plot3(sys_endo_mean_shape(:,1),sys_endo_mean_shape(:,2),sys_endo_mean_shape(:,3),'g.')
plot3(sys_endo_new_shape_max(:,1),sys_endo_new_shape_max(:,2),sys_endo_new_shape_max(:,3),'ro')
plot3(sys_endo_new_shape_min(:,1),sys_endo_new_shape_min(:,2),sys_endo_new_shape_min(:,3),'bo')
legend 'mean' 'max b' 'min b'

figure
title 'PCA - eigenmode visualisation - systolic myocardium'
hold on
plot3(sys_myo_mean(:,1),sys_myo_mean(:,2),sys_myo_mean(:,3),'g.')
plot3(sys_myo_new_shape_2_max(:,1),sys_myo_new_shape_2_max(:,2),sys_myo_new_shape_2_max(:,3),'ro')
plot3(sys_myo_new_shape_2_min(:,1),sys_myo_new_shape_2_min(:,2),sys_myo_new_shape_2_min(:,3),'bo')
legend 'mean' 'max b' 'min b'

figure
for c = -1:0.1:1
plot3D(sys_myo_mean)
hold on
 c   

sys_myo_new_shape(:,1) = sys_myo_mean(:) + sys_myo_principle_eigenvectors(:,1)*c*sys_myo_max_b(1,1);
sys_myo_new_shape(:,2) = sys_myo_mean(:) + sys_myo_principle_eigenvectors(:,2)*c*sys_myo_max_b(2,2);
% 
 plot3D(reshape(sys_myo_new_shape(:,1), [2178 3]))
plot3D(reshape(sys_myo_new_shape(:,2), [2178 3]))


pause
axis equal
hold off
end


%% visualise eigenmodes of shape variation - meshes
figure
patch('Vertices',sys_endo_new_shape_max,'Faces',data(1).systolic.endo.tri,'FaceColor','red')
figure
patch('Vertices',sys_endo_new_shape_min,'Faces',data(1).systolic.endo.tri,'FaceColor','blue')

dia_endo_new_shape_min = reshape(dia_endo_mean_shape + dia_endo_principle_eigenvector.*dia_endo_min_b, [1089 3]);
dia_endo_new_shape_max = reshape(dia_endo_mean_shape + dia_endo_principle_eigenvector.*dia_endo_max_b, [1089 3]);
dia_endo_mean_shape = reshape(dia_endo_mean_shape, [1089 3]);
figure 
hold on
plot3(dia_endo_mean_shape(:,1),dia_endo_mean_shape(:,2),dia_endo_mean_shape(:,3),'g.')
plot3(dia_endo_new_shape_max(:,1),dia_endo_new_shape_max(:,2),dia_endo_new_shape_max(:,3),'ro')
plot3(dia_endo_new_shape_min(:,1),dia_endo_new_shape_min(:,2),dia_endo_new_shape_min(:,3),'bo')
title 'dia endo eigenmode variation'
legend 'mean' 'max b' 'min b'

figure
patch('Vertices',dia_endo_new_shape_max,'Faces',data(1).diastolic.endo.tri,'FaceColor','red')
figure
patch('Vertices',dia_endo_new_shape_min,'Faces',data(1).diastolic.endo.tri,'FaceColor','blue')
 
% ICA
%[icaOut] = fastica(dia_endo_covariance_matrix)


disp('finished shape modeling')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SUPPORT VECTOR MACHINES - PDM 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%7
% ??? should I classify in terms of eigenvalue of each case (found by doing PCA on each case) or by b values,
%    calculated by b = (shape - meanshape)*eigenvectors' ??????????????????
determine_sys_myo_sum_of_shapes = zeros(size(data(2).systolic.myo.xyz(:),1),1);
for i = data(1).DETERMINE_indices'
    for d = find(data(1).DETERMINE_indices==i)
        
    sys_myo_tmp = data(i).systolic.myo.xyz(:);
    determine_sys_myo_sum_of_shapes = sys_myo_tmp + determine_sys_myo_sum_of_shapes;
    determine_sys_myo_mean_shape =  determine_sys_myo_sum_of_shapes/400;
    end
end
% find b values

%%
for i = data(1).DETERMINE_indices'
    for d = find(data(1).DETERMINE_indices==i)
       
%     DETERMINE_dia_endo_b(d,1) = (dia_endo_principle_eigenvectors')*(data(i).diastolic.endo.xyz(:) - dia_endo_mean_shape(:));
%     DETERMINE_sys_endo_b1(d,1) = (sys_endo_principle_eigenvectors(:,1)')*(data(i).systolic.endo.xyz(:) - determine_sys_endo_mean_shape(:));
%     DETERMINE_sys_endo_b2(d,1) = (sys_endo_principle_eigenvectors(:,2)')*(data(i).systolic.endo.xyz(:) - determine_sys_endo_mean_shape(:));
%     DETERMINE_sys_endo_b3(d,1) = (sys_endo_principle_eigenvectors(:,3)')*(data(i).systolic.endo.xyz(:) - determine_sys_endo_mean_shape(:));
%     DETERMINE_sys_endo_b4(d,1) = (sys_endo_principle_eigenvectors(:,4)')*(data(i).systolic.endo.xyz(:) - determine_sys_endo_mean_shape(:));
%     DETERMINE_sys_endo_b5(d,1) = (sys_endo_principle_eigenvectors(:,5)')*(data(i).systolic.endo.xyz(:) - determine_sys_endo_mean_shape(:));
%    
    DETERMINE_sys_myo_b1(d,1) = (sys_myo_principle_eigenvectors(:,1)')*(data(i).systolic.myo.xyz(:) - sys_myo_mean(:));
    DETERMINE_sys_myo_b2(d,1) = (sys_myo_principle_eigenvectors(:,2)')*(data(i).systolic.myo.xyz(:) - sys_myo_mean(:));
    DETERMINE_sys_myo_b3(d,1) = (sys_myo_principle_eigenvectors(:,3)')*(data(i).systolic.myo.xyz(:) - sys_myo_mean(:));
    DETERMINE_sys_myo_b4(d,1) = (sys_myo_principle_eigenvectors(:,4)')*(data(i).systolic.myo.xyz(:) - sys_myo_mean(:));
%     DETERMINE_sys_myo_b5(d,1) = (sys_myo_principle_eigenvectors(:,5)')*(data(i).systolic.myo.xyz(:) - sys_myo_mean(:));
    end
end
% MESA_indices = [MESA_indices(1);MESA_indices(3:end)];
for i = data(1).MESA_indices'
    for d = find(data(1).MESA_indices==i)
        d
%     MESA_dia_endo_b(d,1) = (dia_endo_principle_eigenvectors')*(data(i).diastolic.endo.xyz(:) - dia_endo_mean_shape(:));
%     MESA_sys_endo_b1(d,1) = (sys_endo_principle_eigenvectors(:,1)')*(data(i).systolic.endo.xyz(:) - sys_endo_mean_shape(:));
%     MESA_sys_endo_b2(d,1) = (sys_endo_principle_eigenvectors(:,2)')*(data(i).systolic.endo.xyz(:) - sys_endo_mean_shape(:));
%     MESA_sys_endo_b3(d,1) = (sys_endo_principle_eigenvectors(:,3)')*(data(i).systolic.endo.xyz(:) - sys_endo_mean_shape(:));
%     MESA_sys_endo_b4(d,1) = (sys_endo_principle_eigenvectors(:,4)')*(data(i).systolic.endo.xyz(:) - sys_endo_mean_shape(:));
%     MESA_sys_endo_b5(d,1) = (sys_endo_principle_eigenvectors(:,5)')*(data(i).systolic.endo.xyz(:) - sys_endo_mean_shape(:));

    
    MESA_sys_myo_b1(d,1) = (sys_myo_principle_eigenvectors(:,1)')*(data(i).systolic.myo.xyz(:) - sys_myo_mean(:));
    MESA_sys_myo_b2(d,1) = (sys_myo_principle_eigenvectors(:,2)')*(data(i).systolic.myo.xyz(:) - sys_myo_mean(:));
    MESA_sys_myo_b3(d,1) = (sys_myo_principle_eigenvectors(:,3)')*(data(i).systolic.myo.xyz(:) - sys_myo_mean(:));
    MESA_sys_myo_b4(d,1) = (sys_myo_principle_eigenvectors(:,4)')*(data(i).systolic.myo.xyz(:) - sys_myo_mean(:));
%     MESA_sys_myo_b5(d,1) = (sys_myo_principle_eigenvectors(:,5)')*(data(i).systolic.myo.xyz(:) - sys_myo_mean(:));
    
    end
end

%%
% trainingdata(:,1) = [DETERMINE_sys_myo_b2 ; MESA_sys_myo_b2]; 
% trainingdata(:,2) = [DETERMINE_sys_myo_b3 ; MESA_sys_myo_b3];
% trainingdata(:,3) = [DETERMINE_sys_myo_b5 ; MESA_sys_myo_b5];

trainingdata(:,1) = [DETERMINE_sys_myo_b1 ; MESA_sys_myo_b1]; 
trainingdata(:,2) = [DETERMINE_sys_myo_b2 ; MESA_sys_myo_b2];
trainingdata(:,3) = [DETERMINE_sys_myo_b3 ; MESA_sys_myo_b3];
trainingdata(:,4) = [DETERMINE_sys_myo_b4 ; MESA_sys_myo_b4];
% trainingdata(:,5) = [DETERMINE_sys_myo_b5 ; MESA_sys_myo_b5];
% trainingdata(:,1) = [DETERMINEscore(:,1)  ; MESAscore(:,1)]; 
% trainingdata(:,2) = [DETERMINEscore(:,2) ; MESAscore(:,2)];
% MESA_EF = [MESA_EF(1) ; MESA_EF(3:end)];
% trainingdata(:,3) = [DETERMINEscore(:,3) ; MESAscore(:,3)];
names = char(200,1);
% names(1:100,1) = 'd'; 
% names(101:199,1) = 'm';
names(1:100,1) = 1; 
names(101:200,1) = 2;

% % svmtrain will be removed from later versions of matlab.
% svmStruct = svmtrain(trainingdata,names,'ShowPlot', true, 'kernel_function', 'rbf');
%%
% SVMModel = fitcsvm(trainingdata, group, 'Standardize', true)
SVMModel = fitcsvm(trainingdata(:,1:3), names, 'KernelFunction', 'linear' );
% cross-validation of the SVM
CVSVMModel = crossval(SVMModel);
misclassification_rate = kfoldLoss(CVSVMModel);
classification_rate = 1 - misclassification_rate
%% kmeans
idx = kmeans(trainingdata(:,1:3),2);
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
%%
figure
subplot 231
title 'sys myo b1'
nbins = 30;
hold on
histogram(DETERMINE_sys_myo_b1, nbins)
histogram(MESA_sys_myo_b1, nbins)
hold off


subplot 232
title 'b2'
nbins = 30;
hold on
histogram(DETERMINE_sys_myo_b2, nbins)
histogram(MESA_sys_myo_b2, nbins)
hold off


subplot 233
title 'b3'
nbins = 30;
hold on
histogram(DETERMINE_sys_myo_b3, nbins)
histogram(MESA_sys_myo_b3, nbins)
hold off


subplot 234
title 'b4'
nbins = 30;
hold on
histogram(DETERMINE_sys_myo_b4, nbins)
histogram(MESA_sys_myo_b4, nbins)
hold off


subplot 235
title 'b5'
nbins = 30;
hold on
histogram(DETERMINE_sys_myo_b5, nbins)
histogram(MESA_sys_myo_b5, nbins)
hold off

%%
hold on
xlabel 'b1'
ylabel 'b2'
plot(DETERMINE_sys_myo_b1, DETERMINE_sys_myo_b2, 'o')
plot(MESA_sys_myo_b1, MESA_sys_myo_b2, 'o')
hold off


hold on
xlabel 'b1'
ylabel 'b3'
plot(DETERMINE_sys_myo_b1, DETERMINE_sys_myo_b3, 'o')
plot(MESA_sys_myo_b1, MESA_sys_myo_b3, 'o')
hold off

figure
hold on
xlabel 'b1'
ylabel 'b4'
plot(DETERMINE_sys_myo_b1, DETERMINE_sys_myo_b4, 'o')
plot(MESA_sys_myo_b1, MESA_sys_myo_b4, 'o')
hold off

figure
hold on
xlabel 'b1'
ylabel 'b5'
plot(DETERMINE_sys_myo_b1, DETERMINE_sys_myo_b5, 'o')
plot(MESA_sys_myo_b1, MESA_sys_myo_b5, 'o')
hold off

figure
hold on
xlabel 'b2'
ylabel 'b3'
plot(DETERMINE_sys_myo_b2, DETERMINE_sys_myo_b3, 'o')
plot(MESA_sys_myo_b2, MESA_sys_myo_b3, 'o')
hold off

figure
hold on
xlabel 'b2'
ylabel 'b4'
plot(DETERMINE_sys_myo_b2, DETERMINE_sys_myo_b4, 'o')
plot(MESA_sys_myo_b2, MESA_sys_myo_b4, 'o')
hold off

figure
hold on
xlabel 'b2'
ylabel 'b4'
plot(DETERMINE_sys_myo_b2, DETERMINE_sys_myo_b4, 'o')
plot(MESA_sys_myo_b2, MESA_sys_myo_b4, 'o')
hold off
 %%
figure
hold on
xlabel 'b2'
ylabel 'b3'
zlabel 'b5'
plot3(DETERMINE_sys_myo_b2 , DETERMINE_sys_myo_b3, DETERMINE_sys_myo_b5, 'o')
plot3(MESA_sys_myo_b2, MESA_sys_myo_b3, MESA_sys_myo_b5, 'o')
hold off

figure
title ' sys myo '
hold on
xlabel 'b1'
ylabel 'b2'
zlabel 'b3'
plot3(DETERMINE_sys_myo_b1 , DETERMINE_sys_myo_b2, DETERMINE_sys_myo_b3, 'o')
plot3(MESA_sys_myo_b1, MESA_sys_myo_b2, MESA_sys_myo_b3, 'o')
hold off


classOrder = SVMModel.ClassNames
figure
title 'SVM'
gscatter(trainingdata(:,1),trainingdata(:,2),names)
hold on
sv = SVMModel.SupportVectors;
plot(sv(:,1),sv(:,2),'ko','MarkerSize',10)
legend('DETERMINE','MESA','Support Vector')
% x = 10:80;
% plot(10:80, (0.005*x)+ -7.0407);
hold off

%% PCA using built-in matlab function
% for i = 1
% [coeff,score,latent,tsquared,explained,mu] = pca(data(2).systolic.endo.xyz);
% pcacoeff(:,i) = coeff(:);
% pcascores(:,i) = scores(:);
% end
% 
% MESA_indices = [MESA_indices(1);MESA_indices(3:end)];
% for i = MESA_indices'
%     for d = find(MESA_indices==i)
%     
%     MESAdata(d,1:3267) = data(i).systolic.endo.xyz(:);
%     end
% end
% 
% for i = DETERMINE_indices'
%     for d = find(DETERMINE_indices==i)
%     
%     DETERMINEdata(d,1:3267) = data(i).systolic.endo.xyz(:);
%     end
% end
% combined = [MESA_indices;DETERMINE_indices]'
% for i = combined
%     for d = find(combined==i)
%     testdata(d,1:3267) = data(i).systolic.endo.xyz(:);
%     end
% end
% % [MESAcoeff,MESAscore,MESAlatent,~,MESAexplained]  = pca(MESAdata);
% % [DETERMINEcoeff,DETERMINEscore,DETERMINElatent,~,DETERMINEexplained]  = pca(DETERMINEdata);
% 
% [coeff,score,latent,~,explained]  = pca(testdata);
% hold on
% plot(MESAscore(:,1),MESAscore(:,2),'+')
% plot(DETERMINEscore(:,1),DETERMINEscore(:,2),'o')
% % plot3(score(:,1),score(:,2),score(:,3),'+')
