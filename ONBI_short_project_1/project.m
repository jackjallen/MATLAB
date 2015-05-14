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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATE MYOCARDIUM VOLUMES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% when run after calcVolumes() etc, this gives warnings...


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

plotROC(sensitivity, specificity)
plotAccuracy(accuracy)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATE MESH AREAS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate triangle side lengths
%!!!!! NEED TO CORRECT THIS TO INCLUDE THE LIDS
disp('started calculating triangular mesh areas and triangle side lengths')
[data] = calcTriMeshAreas(data);
disp('finished calculating triangular mesh areas and triangle side lengths')
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
% %initialise the sums of the shapes as zero
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
procrustes_case = i
%     for d = 1:find(concatIndices'==i)
% % transform endo and epi (individually)
% [data(i).diastolic.endo.procrustes_d, data(i).diastolic.endo.xyz] = procrustes(dia_endo_reference, data(i).diastolic.endo.xyz(:));
% [data(i).diastolic.epi.procrustes_d, data(i).diastolic.epi.xyz] = procrustes(dia_epi_reference, data(i).diastolic.epi.xyz(:));
% [data(i).systolic.endo.procrustes_d, data(i).systolic.endo.xyz] = procrustes(sys_endo_reference, data(i).systolic.endo.xyz(:));
% [data(i).systolic.epi.procrustes_d, data(i).systolic.epi.xyz] = procrustes(sys_epi_reference, data(i).systolic.epi.xyz(:));

% transform endo and epi as one shape vector
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SHAPE MODELING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('start shape modeling')
%% PCA
% see Cootes tutorial for help: http://personalpages.manchester.ac.uk/staff/timothy.f.cootes/Models/app_models.pdf

% do PCA on MESA and DETERMINE cases only
for i = 1:400
    for d = 1:size(data(1).DETERMINE_indices,1)
        if data(1).DETERMINE_indices(d)==i
            DETERMINE_diastolic_myo_shapes(d,:) = all_diastolic_myo_shapes(i,:);
            DETERMINE_systolic_myo_shapes(d,:) = all_systolic_myo_shapes(i,:);
        end
        if data(1).MESA_indices(d)==i
            MESA_diastolic_myo_shapes(d,:) = all_diastolic_myo_shapes(i,:);
            MESA_systolic_myo_shapes(d,:) = all_systolic_myo_shapes(i,:);
        
        end
    end
end

training_diastolic_myo_shapes = [DETERMINE_diastolic_myo_shapes ; MESA_diastolic_myo_shapes];
training_systolic_myo_shapes = [DETERMINE_systolic_myo_shapes ; MESA_systolic_myo_shapes];

% find mean shape and then use it to find the covariance matrices
disp('started calculating covariance and mean shape')

dia_myo_cov_mat = cov(training_diastolic_myo_shapes);
sys_myo_cov_mat = cov(training_systolic_myo_shapes);
 
disp('finished calculating covariance and mean shape')
% find eigenvectors of covariance matrix (of just MESA)
disp('started calculating eigenvectors of covariance matrix')

[dia_myo_eigenvectors,dia_myo_eigenvalues]=eig(dia_myo_cov_mat);
[sys_myo_eigenvectors,sys_myo_eigenvalues]=eig(sys_myo_cov_mat);

% % Each eigenvalue gives the variance of the data about the mean in the
% % direction of the corresponding eigenvector. Compute the total variance
%  totalVariance = sum(dia_endo_eigenvalues(:));
% proportion = 1;
% % %Choose the first t largest eigenvalues such that
%  p = proportion*totalVariance;
% a = dia_endo_eigenvalues/sum(dia_endo_eigenvalues(:));
disp('finished calculating eigenvectors of covariance matrix')
%
nEigenvalues = 5;
sorted_dia_myo_eigenvalues = sort(dia_myo_eigenvalues(:), 'descend');
sorted_sys_myo_eigenvalues = sort(sys_myo_eigenvalues(:), 'descend');
principle_dia_myo_eigenvalues = sorted_dia_myo_eigenvalues(1:nEigenvalues);
principle_sys_myo_eigenvalues = sorted_sys_myo_eigenvalues(1:nEigenvalues);

for n = 1:nEigenvalues 
[dia_myo_eRows(1,n), dia_myo_eCols(1,n)] = find(dia_myo_eigenvalues == principle_dia_myo_eigenvalues(n,1));
[sys_myo_eRows(1,n), sys_myo_eCols(1,n)] = find(sys_myo_eigenvalues == principle_sys_myo_eigenvalues(n,1)); 

principle_dia_myo_eigenvectors(:,n) = dia_myo_eigenvectors(:, dia_myo_eCols(1,n));
principle_sys_myo_eigenvectors(:,n) = sys_myo_eigenvectors(:, sys_myo_eCols(1,n));

end 

% Find b values - using +/- 3*sqrt(eigenvalue) as b
% systolic, endocardium
dia_myo_min_b = - 3*sqrt(principle_dia_myo_eigenvalues);
dia_myo_max_b = 3*sqrt(principle_dia_myo_eigenvalues);
sys_myo_min_b = - 3*sqrt(principle_sys_myo_eigenvalues);
sys_myo_max_b = 3*sqrt(principle_sys_myo_eigenvalues);

% Calculate new shapes for visualisation
% mean shape + the sum of each eigenvector times it's value of b
%set of the principle eigenvectors to use (eigenvalue_range)
%(can equal the number_of_eigenvectors already extracted, or we can vary it for plotting)
nEigenmodes = 5;
for n = 1 : nEigenmodes
% dia_endo_new_shape_min(:,n) = dia_endo_mean_shape + dia_endo_principle_eigenvectors(:,eigenvalues(1,n)).*dia_endo_min_b(1,eigenvalues(1,n));
% dia_endo_new_shape_max(:,n) = dia_endo_mean_shape + dia_endo_principle_eigenvectors(:,eigenvalues(1,n)).*dia_endo_max_b(1,eigenvalues(1,n));
% sys_endo_new_shape_min(:,n) = sys_endo_mean_shape + sys_endo_principle_eigenvectors(:,eigenvalues(1,n)).*sys_endo_min_b(1,eigenvalues(1,n));
% sys_endo_new_shape_max(:,n) = sys_endo_mean_shape + sys_endo_principle_eigenvectors(:,eigenvalues(1,n)).*sys_endo_max_b(1,eigenvalues(1,n));
sys_myo_new_shape_min(:,n) = sys_myo_mean(:) + principle_sys_myo_eigenvectors(:,n)*sys_myo_min_b(n,1);
sys_myo_new_shape_max(:,n) = sys_myo_mean(:) + principle_sys_myo_eigenvectors(:,n)*sys_myo_max_b(n,1);
dia_myo_new_shape_min(:,n) = dia_myo_mean(:) + principle_dia_myo_eigenvectors(:,n)*dia_myo_min_b(n,1);
dia_myo_new_shape_max(:,n) = dia_myo_mean(:) + principle_dia_myo_eigenvectors(:,n)*dia_myo_max_b(n,1);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Point distribution model (PDM)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%7

% find mean shapes
mean_diastolic_myo_shape = mean(training_diastolic_myo_shapes);
mean_systolic_myo_shape = mean(training_systolic_myo_shapes);

% find b values for the DETERMINE and MESA cases
for b = 1:5
    for i = data(1).DETERMINE_indices'
        for d = find(data(1).DETERMINE_indices==i)
            
            DETERMINE_dia_myo_b(d,b) = (principle_dia_myo_eigenvectors(:,b)')*(data(i).diastolic.myo.xyz(:) - mean_diastolic_myo_shape');
            DETERMINE_sys_myo_b(d,b) = (principle_sys_myo_eigenvectors(:,b)')*(data(i).systolic.myo.xyz(:) - mean_systolic_myo_shape');
        end
    end
    
    for i = data(1).MESA_indices'
        for d = find(data(1).MESA_indices==i)
            %         d
            MESA_dia_myo_b(d,b) = (principle_dia_myo_eigenvectors(:,b)')*(data(i).diastolic.myo.xyz(:) - mean_diastolic_myo_shape');
            MESA_sys_myo_b(d,b) = (principle_sys_myo_eigenvectors(:,b)')*(data(i).systolic.myo.xyz(:) - mean_systolic_myo_shape');
            
        end
    end
    
end

%% visualise b values - histograms
b = 1
figure
nbins = 30;
hold on
histogram(DETERMINE_sys_myo_b(:,b), nbins)
histogram(MESA_sys_myo_b(:,b), nbins)
hold off

% visualise b values - 2D plot
b = 1:2;
figure
hold on
xlabel 'b2'
ylabel 'b4'
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
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CLASSIFICATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SVM
% set training data
trainingdata(:,1:5) = [DETERMINE_sys_myo_b(:,1:5) ; MESA_sys_myo_b(:,1:5)]; %from the PDM
trainingdata(:,6) = [DETERMINE_EF ; MESA_EF];
names = char(200,1);
% names(1:100,1) = 'd'; 
% names(101:199,1) = 'm';
names(1:100,1) = 1; 
names(101:200,1) = 2;

% svmtrain (will be removed from later versions of matlab)
svmStruct = svmtrain(trainingdata(:,1:2),names,'ShowPlot', true, 'kernel_function', 'linear');

% SVMModel = fitcsvm(trainingdata, group, 'Standardize', true)
SVMModel = fitcsvm(trainingdata(:,:), names, 'KernelFunction', 'linear' );
% cross-validation of the SVM
CVSVMModel = crossval(SVMModel);
misclassification_rate = kfoldLoss(CVSVMModel);
classification_rate = 1 - misclassification_rate

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