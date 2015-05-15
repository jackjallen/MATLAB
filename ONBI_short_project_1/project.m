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

[data, DETERMINE_strokeVolumes, DETERMINE_ejectionFractions,  MESA_strokeVolumes,  MESA_ejectionFractions] = calcEjectionFraction(data);

plotEF(data)
disp('finished calculating ejection fractions')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  CALCULATE ACCURACIES 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('started calculating accuracies')

[data, accuracy, sensitivity, specificity] = calcAccuracyEF(data);

plotROC(sensitivity, specificity)
plotAccuracy(data, accuracy)


disp('finished calculating accuracies')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATE MESH AREAS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('started calculating triangular mesh areas and triangle side lengths')
[data] = calcTriMeshAreas(data);
disp('finished calculating triangular mesh areas and triangle side lengths')

% calculate area to volume ratios
[data] = calcAVR(data);

figure
hold on
histogram(cell2mat({data(data(1).MESA_indices).MESA_sys_myo_AVratio}))
histogram(cell2mat({data(data(1).DETERMINE_indices).DETERMINE_sys_myo_AVratio}))

figure
title 'systolic endocardium, surface area to volume ratio'
hold on
histogram(cell2mat({data(data(1).MESA_indices).MESA_sys_endo_AVratio}))
histogram(cell2mat({data(data(1).DETERMINE_indices).DETERMINE_sys_endo_AVratio}))

[data, accuracy, sensitivity, specificity] = calcAccuracyAVratio(data);
plotROC(sensitivity, specificity)
plotAccuracy(data, accuracy)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GENERAL PROCRUSTES ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is done by concatenating the endo and epi shape vectors to produce
% longer vectors (myo.xyz). Procrustes analysis is then performed on the
% two myo.xyz vectors. Finally, the transformed myo.xyz vectors are split
% to extract transformed endo and epi shape vectors.
disp('starting procrustes analysis')

%store dimensions of the shapes
[shape_nRows , shape_nCols] = size(data(2).diastolic.endo.xyz(1:1089,:));

%use the initial data clouds as references
initial_reference_id = 401; %sets which case will be the initial reference

%remove the B points (added in calcVolumes())
dia_endo_reference = data(initial_reference_id).diastolic.endo.xyz(1:1089,:);
dia_epi_reference = data(initial_reference_id).diastolic.epi.xyz(1:1089,:);
sys_endo_reference = data(initial_reference_id).systolic.endo.xyz(1:1089,:);
sys_epi_reference = data(initial_reference_id).systolic.epi.xyz(1:1089,:);

dia_myo_reference = reshape([dia_endo_reference(:) ; dia_epi_reference(:)], [2*1089 3]);
sys_myo_reference = reshape([sys_endo_reference(:)  ; sys_epi_reference(:)], [2*1089 3]);

% Think of endo and epi as one shape (concatenate).
for i = 1:401 %SMM001 has already been replaced by SMM401 in the script
data(i).diastolic.myo.xyz = [data(i).diastolic.endo.xyz(1:1089,:) ; data(i).diastolic.epi.xyz(1:1089,:)];
data(i).systolic.myo.xyz = [data(i).systolic.endo.xyz(1:1089,:) ; data(i).systolic.epi.xyz(1:1089,:)];
end

%p = number of times procrustes is performed.
procrustes_iterations = 3;
[data, dia_myo_mean, sys_myo_mean, MESA_diastolic_myo_shapes, DETERMINE_systolic_myo_shapes, all_training_diastolic_myo_shapes, all_training_systolic_myo_shapes ] = calcProcrustes(data, procrustes_iterations, dia_myo_reference, sys_myo_reference);

% % reshape individual endo and epi shapes.
% dia_endo_mean = reshape(dia_endo_mean,size(data(1).diastolic.endo.xyz));
% dia_epi_mean = reshape(dia_epi_mean,size(data(1).diastolic.endo.xyz));
% sys_endo_mean = reshape(sys_endo_mean,size(data(1).diastolic.endo.xyz));
% sys_epi_mean = reshape(sys_epi_mean,size(data(1).diastolic.endo.xyz));
dia_myo_reference = reshape(dia_myo_mean, [2*shape_nRows, shape_nCols]);
sys_myo_reference = reshape(sys_myo_mean, [2*shape_nRows, shape_nCols]);

[data] = extractEndoEpi(data, shape_nRows);

disp('finished procrustes analysis')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SHAPE MODELING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PCA
% see Cootes tutorial for help: http://personalpages.manchester.ac.uk/staff/timothy.f.cootes/Models/app_models.pdf
disp('started PCA')
% do PCA on MESA and DETERMINE cases only
% nCases = 401;
% for i = 1:nCases 
for d = data(1).DETERMINE_indices'
 
    DETERMINE_diastolic_myo_shapes(d,:) = all_training_diastolic_myo_shapes(d,:);
    DETERMINE_systolic_myo_shapes(d,:) = all_training_systolic_myo_shapes(d,:);
   
end
for m = data(1).MESA_indices'
    MESA_diastolic_myo_shapes(m,:) = all_training_diastolic_myo_shapes(m,:);
    MESA_systolic_myo_shapes(m,:) = all_training_systolic_myo_shapes(m,:);
    
end
% delete empty rows
DETERMINE_diastolic_myo_shapes( ~any(DETERMINE_diastolic_myo_shapes,2), : ) = [];
DETERMINE_systolic_myo_shapes( ~any(DETERMINE_systolic_myo_shapes,2), : ) = [];
MESA_diastolic_myo_shapes( ~any(MESA_diastolic_myo_shapes,2), : ) = [];
MESA_systolic_myo_shapes( ~any(MESA_systolic_myo_shapes,2), : ) = [];

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

disp('finished PCA')
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

%% visualise modes

f1= figure;
% screen_size = get(0, 'ScreenSize');
% set(f1,'Position', [0 0 screen_size(3) screen_size(4)])
pause
for n = 1:2
    for c = -1:0.1:1
        
        sys_myo_new_shape(:,1) = sys_myo_mean(:) + principle_sys_myo_eigenvectors(:,1)*c*sys_myo_max_b(1,1);
        sys_myo_new_shape(:,2) = sys_myo_mean(:) + principle_sys_myo_eigenvectors(:,2)*c*sys_myo_max_b(2,1);
        %
        
        subplot 121
        plot3D(reshape(sys_myo_mean,[2178 , 3]))
        set(gca ,'XLim',[-50, 10], 'YLim', [-50 10], 'ZLim', [-80, 0])
        xlabel 'x', ylabel 'y', zlabel 'z'
        hold on
        plot3D(reshape(sys_myo_new_shape(:,1), [2178 3]))
        set(gca ,'XLim',[-50, 10], 'YLim', [-50 10], 'ZLim', [-80, 30])
        hold off
        legend 'mean' 'mode1'
        
        subplot 122
        plot3D(reshape(sys_myo_mean,[2178 , 3]))
        set(gca ,'XLim',[-50, 10], 'YLim', [-50 10], 'ZLim', [-80, 0])
        xlabel 'x', ylabel 'y', zlabel 'z'
        hold on
        plot3D(reshape(sys_myo_new_shape(:,2), [2178 3]))
        set(gca ,'XLim',[-50, 10], 'YLim', [-50 10], 'ZLim', [-80,30])
        legend 'mean' 'mode2'
        pause(0.3)
        axis equal
        
        hold off
    end
end
%% visualise b values - histograms
b = 5;
figure
nbins = 30;
hold on
histogram(DETERMINE_sys_myo_b(:,b), nbins)
histogram(MESA_sys_myo_b(:,b), nbins)
hold off

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
 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CLASSIFICATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SVM
% set training data
trainingdata(:,1:5) = [DETERMINE_sys_myo_b(:,1:5) ; MESA_sys_myo_b(:,1:5)]; %from the PDM
trainingdata(:,6) = [DETERMINE_ejectionFractions' ; MESA_ejectionFractions'];
names = char(200,1);
% names(1:100,1) = 'd'; 
% names(101:199,1) = 'm';
names(1:100,1) = 1; 
names(101:200,1) = 2;

b = [1;4;6]
figure
title 'SVM classification'
hold on 
plot3D(trainingdata(1:100,b))
plot3D(trainingdata(101:200,b))
hold off
legend ' DETERMINE' 'MESA'
xlabel 'b1'
ylabel 'b4'
zlabel 'Ejection fractions'

% svmtrain (will be removed from later versions of matlab)
svmStruct = svmtrain(trainingdata(:,[4;6]),names,'ShowPlot', true, 'kernel_function', 'linear');
% 

% %  which is the best trio of b values for classification?
% tic
% for param1 = 1:5
%     for param2 = 1:5
%         for param3 = 1:5
%             % SVMModel = fitcsvm(trainingdata, group, 'Standardize', true)
%             SVMModel = fitcsvm(trainingdata(:,[param1;param2;param3]), names, 'KernelFunction', 'linear' );
%             % cross-validation of the SVM
%             CVSVMModel = crossval(SVMModel);
%             misclassification_rate = kfoldLoss(CVSVMModel);
%             classification_rates(param1,param2,param3) = 1 - misclassification_rate
%         end
%     end
% end
% toc

tic
for param1 = 1:5
    for param2 = 1:5
        param1
        param2
        % SVMModel = fitcsvm(trainingdata, group, 'Standardize', true)
        SVMModel = fitcsvm(trainingdata(:,[param1;param2]), names, 'KernelFunction', 'linear' );
        % cross-validation of the SVM
        CVSVMModel = crossval(SVMModel);
        misclassification_rate = kfoldLoss(CVSVMModel);
        classification_rates(param1,param2) = 1 - misclassification_rate                
    end 
end
toc
tic
for param1 = 1:5
    for param2 = 1:5
        param1
        param2
        % SVMModel = fitcsvm(trainingdata, group, 'Standardize', true)
        SVMModel = fitcsvm(trainingdata(:,[param1;param2;6]), names, 'KernelFunction', 'linear' );
        % cross-validation of the SVM
        CVSVMModel = crossval(SVMModel);
        misclassification_rate = kfoldLoss(CVSVMModel);
        classification_rates_withEF(param1,param2) = 1 - misclassification_rate                
    end 
end
toc
imagesc(classification_rates)
imagesc(classification_rates_withEF)

figure
subplot 121
% surf(classification_rates)
imagesc(classification_rates)
title ' classification rates'
% xlim ([1 ; 5]')
% ylim ([1 ; 5])
xlabel 'b'
ylabel 'b'
caxis([.5, .9])
colorbar

subplot 122
% surf(classification_rates_withEF(:,:,6))
imagesc(classification_rates_withEF(:,:,6))
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