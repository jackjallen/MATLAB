%% ONBI Short project 1
% Jack Allen
% Supervisor: Vicente Grau
%
clear all
close all
clc

setenv('path',[getenv('path'),';','C:\Users\jesu2687\Documents\MATLAB\ErnestoCode\Tools\MESHES\vtk_libs']);
addpath('C:\Users\jesu2687\Documents\MATLAB\ONBI_short_project_1')
addpath C:\Users\jesu2687\Documents\MATLAB\ErnestoCode\Tools\MESHES\
addpath C:\Users\jesu2687\Documents\MATLAB\ErnestoCode\Tools\

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
%store volumes
for i = 400
    transformed_data(1).diastolic.endo.volume = dia_endo_volumes(i,1);
    transformed_data(1).diastolic.epi.volume = dia_epi_volumes(i,1);
    transformed_data(1).systolic.endo.volume = sys_endo_volumes(i,1);
    transformed_data(1).systolic.epi.volume = sys_epi_volumes(i,1);
end
%% Myocardium volumes
for i = 1:400 %all patients
 %reshape from a vector to a matrix.
%  %diastolic
 transformed_data(i).diastolic.endo.xyz = reshape(transformed_data(i).diastolic.endo.xyz, [1089, 3]);
 transformed_data(i).diastolic.epi.xyz = reshape(transformed_data(i).diastolic.epi.xyz, [1089, 3]);
 %systolic
 transformed_data(i).systolic.endo.xyz = reshape(transformed_data(i).systolic.endo.xyz, [1089, 3]);
 transformed_data(i).systolic.epi.xyz = reshape(transformed_data(i).systolic.epi.xyz, [1089, 3]);

% vtkCleanPolyData(EPI_ED) fix the possible replicated nodes and spurious
% edges.
% diastolic
transformed_data(i).diastolic.endo.B = vtkFeatureEdges( vtkCleanPolyData(transformed_data(i).diastolic.endo) , 'BoundaryEdgesOn' , [] , 'FeatureEdgesOff' , [] );  %extracting the boundary
transformed_data(i).diastolic.endo.B.xyz = transformed_data(i).diastolic.endo.B.xyz( [2 1 3:end], : );  %fixing the connectivity.
transformed_data(i).diastolic.epi.B = vtkFeatureEdges( vtkCleanPolyData(transformed_data(i).diastolic.epi) , 'BoundaryEdgesOn' , [] , 'FeatureEdgesOff' , [] );  %extracting the boundary
transformed_data(i).diastolic.epi.B.xyz = transformed_data(i).diastolic.epi.B.xyz( [2 1 3:end], : );  %fixing the connectivity.
%systolic
transformed_data(i).systolic.endo.B = vtkFeatureEdges( vtkCleanPolyData(transformed_data(i).systolic.endo) , 'BoundaryEdgesOn' , [] , 'FeatureEdgesOff' , [] );  %extracting the boundary
transformed_data(i).systolic.endo.B.xyz = transformed_data(i).systolic.endo.B.xyz( [2 1 3:end], : );  %fixing the connectivity.
transformed_data(i).systolic.epi.B = vtkFeatureEdges( vtkCleanPolyData(transformed_data(i).systolic.epi) , 'BoundaryEdgesOn' , [] , 'FeatureEdgesOff' , [] );  %extracting the boundary
transformed_data(i).systolic.epi.B.xyz = transformed_data(i).systolic.epi.B.xyz( [2 1 3:end], : );  %fixing the connectivity.


%load a struct containing a manually produced myoB.tri
load('C:\Users\jesu2687\Documents\MATLAB\ONBI_short_project_1\myoB.mat')

%join the list of coordinates for endo and epi to be used to cover the myo.
%diastolic
transformed_data(i).diastolic.myoB.xyz = [transformed_data(i).diastolic.endo.B.xyz ; transformed_data(i).diastolic.epi.B.xyz];
%systolic
transformed_data(i).systolic.myoB.xyz = [transformed_data(i).systolic.endo.B.xyz ; transformed_data(i).systolic.epi.B.xyz];

%append epi surface, endo surface and myo lid
%diastolic
Last_prevID1  = size( data(1).diastolic.endo.xyz , 1 );
Last_prevID2  = size( data(1).diastolic.epi.xyz , 1 );
Last_prevID3 = size(transformed_data(i).diastolic.myoB.xyz,1);
transformed_data(i).diastolic.myoB_full.xyz = [ transformed_data(i).diastolic.endo.xyz ; transformed_data(i).diastolic.myoB.xyz ; transformed_data(i).diastolic.epi.xyz ];
transformed_data(i).diastolic.myoB_full.tri = [ transformed_data(i).diastolic.endo.tri; myoB.tri + Last_prevID1 ; transformed_data(i).diastolic.epi.tri + Last_prevID1 + Last_prevID3 ];
%systolic
Last_prevID1  = size( data(1).systolic.endo.xyz , 1 );
Last_prevID2  = size( data(1).systolic.epi.xyz , 1 );
Last_prevID3 = size(transformed_data(i).systolic.myoB.xyz,1);
transformed_data(i).systolic.myoB_full.xyz = [ transformed_data(i).systolic.endo.xyz ; transformed_data(i).systolic.myoB.xyz ; transformed_data(i).systolic.epi.xyz ];
transformed_data(i).systolic.myoB_full.tri = [ transformed_data(i).systolic.endo.tri; myoB.tri + Last_prevID1 ; transformed_data(i).systolic.epi.tri + Last_prevID1 + Last_prevID3 ];

%make sure that every triangle points outwards.
%diastolic
transformed_data(i).diastolic.myoB_full = FixNormals(transformed_data(i).diastolic.myoB_full);
%systolic
transformed_data(i).systolic.myoB_full = FixNormals(transformed_data(i).systolic.myoB_full);

%size(transformed_data(i).diastolic.myoB_full)

%calculate volume of myocardium.
%diastolic
[transformed_data(i).diastolic.myoVolume, transformed_data(i).diastolic.myoCenterOfMass] = MeshVolume( transformed_data(i).diastolic.myoB_full );
%systolic
[transformed_data(i).systolic.myoVolume, transformed_data(i).systolic.myoCenterOfMass] = MeshVolume( transformed_data(i).systolic.myoB_full );

% %transformed_data(i).diastolic.myodifference_volume = prod(diff( BBMesh( transformed_data(i).diastolic.myoB_full ) , 1  , 1 ) ) - transformed_data(i).diastolic.myoVolume ;   %%it shoud be positive!!
% transformed_data(i).diastolic.myodifference_volume =  transformed_data(i).diastolic.epi.volume - transformed_data(i).diastolic.myoVolume ;   %%it shoud be positive!!
% transformed_data(i).systolic.myodifference_volume = transformed_data(i).systolic.epi.volume - transformed_data(i).systolic.myoVolume ;   %%it shoud be positive!!
end

%% myocardium visualisation
% hold on
% plot3(myoB.xyz(:,1),myoB.xyz(:,2), myoB.xyz(:,3),'+')
% plot3(myoB_full.xyz(:,1),myoB_full.xyz(:,2), myoB_full.xyz(:,3),'o')

% plot3(transformed_data(i).diastolic.myoB.xyz(:,1),transformed_data(i).diastolic.myoB.xyz(:,2), transformed_data(i).diastolic.myoB.xyz(:,3))

% % visualise labelled points on the edge of the myocardium lid
% plot3(transformed_data(i).diastolic.myoB.xyz(:,1) , transformed_data(i).diastolic.myoB.xyz(:,2) ,  transformed_data(i).diastolic.myoB.xyz(:,3) , '.r' )
% text(transformed_data(i).diastolic.myoB.xyz(:,1) ,  transformed_data(i).diastolic.myoB.xyz(:,2) ,  transformed_data(i).diastolic.myoB.xyz(:,3) , arrayfun( @(id)sprintf('%d',id) , 1:size(transformed_data(i).diastolic.myoB.xyz,1) , 'un',false ) )

% %visualise center of mass within myocardium mesh
% cla
% patch('vertices',transformed_data(i).diastolic.myoB_full.xyz,'faces',transformed_data(i).diastolic.myoB_full.tri,'facecolor','g','facealpha',0.1);
% hold on;
% plot3( CenterOfMass(1) , CenterOfMass(2) , CenterOfMass(3) , '*r','markers',20 ); hold off


% % visualise myocardium mesh
% cla
% figure
% patch('vertices',transformed_data(i).diastolic.myoB_full.xyz,'faces',transformed_data(i).diastolic.myoB_full.tri,'facecolor','red')

% %visualise surfaces: endo, epi and myo lid
% cla
% patch('vertices',data(1).diastolic.endo.xyz,'faces',data(1).diastolic.endo.tri,'facecolor','none','EdgeColor','red')
% patch('vertices',data(1).diastolic.epi.xyz,'faces',data(1).diastolic.epi.tri,'facecolor','none','EdgeColor','blue')
% patch('vertices',myoB.xyz,'faces',myoB.tri,'facecolor','g')

%% Shape modeling
[dia_endo_covariance_matrix, dia_endo_mean_shape] = calcCovarianceMatrix(transformed_data);
%[covariance_matrix] = cov(transformed_data);

% PCA
%find eigenvectors for the t largest eigen values.
%store eigenvectors in matrix phi.
%columns of dia_endo_eigenvectors are the eigenvectors
[dia_endo_eigenvectors, dia_endo_eigenvalues] = eig(dia_endo_covariance_matrix);

phi = dia_endo_eigenvectors';
phi = (sortrows(phi))';

% %reverse the eigenvalue matrix
%  dia_endo_eigenvalues = dia_endo_eigenvalues(end:-1:1);
% % reverse the columns
% dia_endo_eigenvectors = dia_endo_eigenvectors(:,end:-1:1); 
%  dia_endo_eigenvectors=dia_endo_eigenvectors';

% model parameters 'b'

%phi = dia_endo_eigenvectors(:);

% select some of the largest eigenvalues
%eigenvalue_indices = find(phi>0.6, 20);
%phi(eigenvalue_indices);

sortedPhi = sort(phi, 'descend');
phi = reshape((sortedPhi.'), size(dia_endo_covariance_matrix)); %each column shows an eigenvetor
%coeff = pca(covariance_matrix)

% calculate new shapes, using different parameters b
%!!!!!!!!!!!WHAT VALUES SHOULD 'b' BE?
% new_shape = mean_shape + (phi)*(b);
for i = 1:400
b(i,1) = (phi.')*(transformed_data(1).diastolic.endo.xyz(:) - dia_endo_mean_shape); %transpose phi
end

% plot in pca space... not sure if this is correct...
for i = 1:400
pc1 = eigenvectors * transformed_data(i).diastolic.endo.xyz;
% pc2 = eigenvectors * transformed_data(i).systolic.endo.xyz;
hold on
plot(pc1(1,:),pc1(2,:),'o');
% plot(pc2(1,:),pc2(2,:),'.');

end

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
