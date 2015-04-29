%% ONBI Short project 1
% Jack Allen
% Supervisor: Vicente Grau
%
clear all
close all
clc
%%
addpath('C:\Users\jesu2687\Documents\MATLAB\ONBI_short_project_1')
% read data (adapted from VG code)
load('C:\Users\jesu2687\Documents\MATLAB\ONBI_short_project_1\project1_data.mat');

%% Procrustes analysis
% subplot 121
% plot3(data(1).diastolic.endo.xyz(:,1),data(1).diastolic.endo.xyz(:,2),data(1).diastolic.endo.xyz(:,3),'+')
% axis equal; axis tight;
% hold on
% plot3(data(2).diastolic.endo.xyz(:,1),data(2).diastolic.endo.xyz(:,2),data(2).diastolic.endo.xyz(:,3),'o')

%use the initial data clouds as references
dia_endo_reference = data(1).diastolic.endo.xyz(:);
dia_epi_reference = data(1).diastolic.epi.xyz(:);
sys_endo_reference = data(1).systolic.endo.xyz(:);
sys_epi_reference = data(1).systolic.epi.xyz(:);

dia_endo_sum = zeros(size(dia_endo_reference));
% for t = 1:2

for i = 1:400
% tmpdata(i).diastolic.endo.xyz = data(i).diastolic.endo.xyz;

   [transformed_data(i).diastolic.endo.procrustes_d, transformed_data(i).diastolic.endo.xyz] = procrustes(dia_endo_reference, data(i).diastolic.endo.xyz(:));
   [transformed_data(i).diastolic.epi.procrustes_d, transformed_data(i).diastolic.epi.xyz] = procrustes(dia_epi_reference, data(i).diastolic.epi.xyz(:));
   [transformed_data(i).systolic.endo.procrustes_d, transformed_data(i).systolic.endo.xyz] = procrustes(sys_endo_reference, data(i).systolic.endo.xyz(:));
   [transformed_data(i).systolic.epi.procrustes_d, transformed_data(i).systolic.epi.xyz] = procrustes(sys_epi_reference, data(i).systolic.epi.xyz(:));

   transformed_data(i).diastolic.endo.tri = data(1).diastolic.endo.tri;
   transformed_data(i).diastolic.epi.tri = data(1).diastolic.epi.tri;
   transformed_data(i).systolic.endo.tri = data(1).systolic.endo.tri;
   transformed_data(i).systolic.epi.tri = data(1).systolic.epi.tri;


 dia_endo_sum = transformed_data(i).diastolic.endo.xyz + dia_endo_sum;
 
end
dia_endo_mean = dia_endo_sum/400;
% sum(dia_endo_reference)
% tmpdata(i).diastolic.endo.xyz  = transformed_data(1).diastolic.endo.xyz;
%end

%% Visualising to check procrustes output
dia_endo_reference = reshape(dia_endo_reference, size(data(1).diastolic.endo.xyz));
transformed_data(2).diastolic.endo.xyz = reshape(transformed_data(2).diastolic.endo.xyz, size(data(1).diastolic.endo.xyz));

hold on
%patient 1
plot3(dia_endo_reference(:,1), dia_endo_reference(:,2), dia_endo_reference(:,3),'o');
%patient 2 before procrustes
plot3(data(2).diastolic.endo.xyz(:,1), data(2).diastolic.endo.xyz(:,2), data(2).diastolic.endo.xyz(:,3),'o'); 
%patient 2 after procrustes
plot3(transformed_data(2).diastolic.endo.xyz(:,1), transformed_data(2).diastolic.endo.xyz(:,2), transformed_data(2).diastolic.endo.xyz(:,3),'o');
legend 'pat1 (reference)' 'pat2 dia endo before procrustes' 'pat2 dia endo after procrustes'

%% make lid for myocardium

%% Calculate volumes
%!!!!!!!!!!!!SPLIT calcVolumes INTO MULTIPLE FUNCTIONS!!!!!!!!!!!
%[transformed_data.systolic.epi.tri, transformed_data.systolic.endo.tri, transformed_data.diastolic.epi.tri, transformed_data.diastolic.endo.tri ] = addLid(transformed_data);

[sys_epi_volumes, sys_endo_volumes, dia_epi_volumes, dia_endo_volumes,data] = calcVolumes(transformed_data);


%% Shape modeling
[covariance_matrix] = calcCovarianceMatrix(transformed_data);
%[covariance_matrix] = cov(transformed_data);

% PCA
%find eigenvectors for the t largest eigen values.
%store eigenvectors in matrix phi.
[eigenvectors, eigenvalues] = eig(covariance_matrix);

plot(eigenvectors(:,1))

eigenvalues = eigenvalues(end:-1:1);
eigenvectors = eigenvectors(:,end:-1:1); 
eigenvectors=eigenvectors';

%plot in pca space... not sure if this is correct...
for i = 1:400
pc1 = eigenvectors * transformed_data(i).diastolic.endo.xyz;
% pc2 = eigenvectors * transformed_data(i).systolic.endo.xyz;
hold on
plot(pc1(1,:),pc1(2,:),'o');
% plot(pc2(1,:),pc2(2,:),'.');

end
phi = eigenvectors(:);
% select some of the largest eigenvalues
%eigenvalue_indices = find(phi>0.6, 20);
%phi(eigenvalue_indices);
sortedPhi = sort(phi, 'descend');
phi = reshape((sortedPhi.'), size(covariance_matrix)); %each column shows an eigenvetor
%coeff = pca(covariance_matrix)
%calculate new shapes, using different parameters b
%!!!!!!!!!!!WHAT VALUES SHOULD 'b' BE?
new_shape = mean_shape + (phi)*(b);
b = (phi.')*(new_shape - mean_shape); %transpose phi

%% calculate triangle side lengths
%!!!!! NEED TO CORRECT THIS TO INCLUDE THE LIDS
endo_sides = calcTriSides(trifac, endo);
epi_sides = calcTriSides(trifac, endo);
%% calculate total areas(heron's forumla)
%!!!!! NEED TO CORRECT THIS TO INCLUDE THE LIDS
endo_area = calcTriMeshArea(endo_sides);
epi_area = calcTriMeshArea(epi_sides);
%% calculate coordinates of the centroids
%!!!!! NEED TO CORRECT THIS...ALREADT DONE BY 'calcVolume'?
endo_centroid = mean(endo);
epi_centroid = mean(epi);
