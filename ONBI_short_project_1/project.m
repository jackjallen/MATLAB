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
addpath C:\Users\jesu2687\Documents\MATLAB\ErnestoCode\Tools
addpath C:\Users\jesu2687\Documents\MATLAB\ErnestoCode\

% read data (adapted from VG code)
load('C:\Users\jesu2687\Documents\MATLAB\ONBI_short_project_1\project1_data.mat');

%% Procrustes analysis
disp('start procrustes analysis')

% subplot 121
% plot3(data(1).diastolic.endo.xyz(:,1),data(1).diastolic.endo.xyz(:,2),data(1).diastolic.endo.xyz(:,3),'+')
% axis equal; axis tight;
% hold on
% plot3(data(2).diastolic.endo.xyz(:,1),data(2).diastolic.endo.xyz(:,2),data(2).diastolic.endo.xyz(:,3),'o')

%use the initial data clouds as references (individual)
dia_endo_reference = data(1).diastolic.endo.xyz(:);
dia_epi_reference = data(1).diastolic.epi.xyz(:);
sys_endo_reference = data(1).systolic.endo.xyz(:);
sys_epi_reference = data(1).systolic.epi.xyz(:);

%make a reference from the endo and epi surfaces (concatenate)
dia_myo_reference = [data(1).diastolic.endo.xyz(:) ; data(1).diastolic.epi.xyz(:)];
sys_myo_reference = [data(1).systolic.endo.xyz(:) ; data(1).systolic.epi.xyz(:)];


dia_endo_sum = zeros(size(dia_endo_reference));
% for t = 1:2
for t = 1:2
for i = 1:400
% tmpdata(i).diastolic.endo.xyz = data(i).diastolic.endo.xyz;

% transform endo and epi, individually
[data(i).diastolic.endo.procrustes_d, data(i).diastolic.endo.xyz] = procrustes(dia_endo_reference, data(i).diastolic.endo.xyz(:));
[data(i).diastolic.epi.procrustes_d, data(i).diastolic.epi.xyz] = procrustes(dia_epi_reference, data(i).diastolic.epi.xyz(:));
[data(i).systolic.endo.procrustes_d, data(i).systolic.endo.xyz] = procrustes(sys_endo_reference, data(i).systolic.endo.xyz(:));
[data(i).systolic.epi.procrustes_d, data(i).systolic.epi.xyz] = procrustes(sys_epi_reference, data(i).systolic.epi.xyz(:));


% transform endo and epi, concatenated

% data(i).diastolic.myo.xyz = [data(i).diastolic.endo.xyz ; data(i).diastolic.epi.xyz];
% data(i).systolic.myo.xyz = [data(i).systolic.endo.xyz ; data(i).systolic.epi.xyz];
% 
[data(i).diastolic.myo.procrustes_d, data(i).diastolic.myo.xyz] = procrustes(dia_myo_reference, data(i).diastolic.myo.xyz(:));
[data(i).systolic.myo.procrustes_d, data(i).systolic.myo.xyz] = procrustes(sys_myo_reference, data(i).systolic.myo.xyz(:));

%
dia_endo_sum = data(i).diastolic.endo.xyz + dia_endo_sum;

% reshape from a vector to a matrix.
%diastolic
data(i).diastolic.endo.xyz = reshape(data(i).diastolic.endo.xyz, size(data(1).diastolic.endo.xyz));
data(i).diastolic.epi.xyz = reshape(data(i).diastolic.epi.xyz, size(data(1).diastolic.endo.xyz));
%systolic
data(i).systolic.endo.xyz = reshape(data(i).systolic.endo.xyz, size(data(1).diastolic.endo.xyz));
data(i).systolic.epi.xyz = reshape(data(i).systolic.epi.xyz, size(data(1).diastolic.endo.xyz));


end
dia_endo_mean = dia_endo_sum/400;
dia_endo_reference = dia_endo_mean;
end
dia_endo_mean = reshape(dia_endo_mean,size(data(1).diastolic.endo.xyz));

disp('finish procrustes analysis')
%% Visualising to check procrustes output
disp('visualise procrustes output')

plotMESH = @(M,varargin)patch('vertices',M.xyz,'faces',M.tri,'edgecolor','k','facecolor','b',varargin{:});

%from vector to matrix
dia_endo_reference = reshape(dia_endo_reference, size(data(1).diastolic.endo.xyz));
data(i).diastolic.endo.xyz = reshape(data(i).diastolic.endo.xyz, size(data(1).diastolic.endo.xyz));

hold on
%patient 1
plot3(dia_endo_reference(:,1), dia_endo_reference(:,2), dia_endo_reference(:,3),'o');
%patient i before procrustes
plot3(data(i).diastolic.endo.xyz(:,1), data(i).diastolic.endo.xyz(:,2), data(i).diastolic.endo.xyz(:,3),'o'); 
%patient i after procrustes
plot3(data(i).diastolic.endo.xyz(:,1), data(i).diastolic.endo.xyz(:,2), data(i).diastolic.endo.xyz(:,3),'o');
%new mean after procrustes
plot3(data(i).diastolic.endo.xyz(:,1), data(i).diastolic.endo.xyz(:,2), data(i).diastolic.endo.xyz(:,3),'+');
%patient i after another procrustes (using new mean as reference)
plot3(data(i).diastolic.endo.xyz(:,1), data(i).diastolic.endo.xyz(:,2), data(i).diastolic.endo.xyz(:,3),'.');

legend 'pat1 (reference)' 'pat2 dia endo before procrustes' 'pat2 dia endo after procrustes' 'new mean after procrustes' 'pat2 dia endo after another procrustes (using new mean as reference)'

%% make lid for myocardium
disp('make lid for myocardium')
%% endo and epi volumes
disp('start calculating endo and epi volumes')
%!!!!!!!!!!!!SPLIT calcVolumes INTO MULTIPLE FUNCTIONS!!!!!!!!!!!
%[transformed_data.systolic.epi.tri, transformed_data.systolic.endo.tri, transformed_data.diastolic.epi.tri, transformed_data.diastolic.endo.tri ] = addLid(transformed_data);

% [sys_epi_volumes, sys_endo_volumes, dia_epi_volumes, dia_endo_volumes] = calcVolumes(transformed_data);
[sys_epi_volumes, sys_endo_volumes, dia_epi_volumes, dia_endo_volumes] = calcVolumes(data);
%store volumes
for i = 400
    data(1).diastolic.endo.volume = dia_endo_volumes(i,1);
    data(1).diastolic.epi.volume = dia_epi_volumes(i,1);
    data(1).systolic.endo.volume = sys_endo_volumes(i,1);
    data(1).systolic.epi.volume = sys_epi_volumes(i,1);
    
    data(1).diastolic.endo.volume = dia_endo_volumes(i,1);
    data(1).diastolic.epi.volume = dia_epi_volumes(i,1);
    data(1).systolic.endo.volume = sys_endo_volumes(i,1);
    data(1).systolic.epi.volume = sys_epi_volumes(i,1);
    
    
    
end
%% plot volume histograms
hold on
nbins = 100;
histogram(dia_endo_volumes,nbins)
histogram(dia_epi_volumes,nbins)
histogram(sys_endo_volumes,nbins)
histogram(sys_epi_volumes,nbins)
legend ' diastolic endo' ' diastolic epi' 'systolic endo' 'systolic epi'
title 'volume histograms'
xlabel 'volume'
ylabel 'frequency'

disp('finish calculating endo and epi volumes')
%% Myocardium volumes
disp('started calculating myocardium volumes')
% for i = 1:400 %all patients

for i = 1:400;

% vtkCleanPolyData(EPI_ED) fix the possible replicated nodes and spurious
% edges.
% diastolic

data(i).diastolic.endo.B = vtkFeatureEdges( vtkCleanPolyData(data(i).diastolic.endo) , 'BoundaryEdgesOn' , [] , 'FeatureEdgesOff' , [] );  %extracting the boundary
data(i).diastolic.endo.B.xyz = data(i).diastolic.endo.B.xyz( [2 1 3:end], : );  %fixing the connectivity.
data(i).diastolic.epi.B = vtkFeatureEdges( vtkCleanPolyData(data(i).diastolic.epi) , 'BoundaryEdgesOn' , [] , 'FeatureEdgesOff' , [] );  %extracting the boundary
data(i).diastolic.epi.B.xyz = data(i).diastolic.epi.B.xyz( [2 1 3:end], : );  %fixing the connectivity.


data(i).systolic.endo.B = vtkFeatureEdges( vtkCleanPolyData(data(i).systolic.endo) , 'BoundaryEdgesOn' , [] , 'FeatureEdgesOff' , [] );  %extracting the boundary
data(i).systolic.endo.B.xyz = data(i).systolic.endo.B.xyz( [2 1 3:end], : );  %fixing the connectivity.
data(i).systolic.epi.B = vtkFeatureEdges( vtkCleanPolyData(data(i).systolic.epi) , 'BoundaryEdgesOn' , [] , 'FeatureEdgesOff' , [] );  %extracting the boundary
data(i).systolic.epi.B.xyz = data(i).systolic.epi.B.xyz( [2 1 3:end], : );  %fixing the connectivity.
% find nearest points on endo and epi boundaries and assign them as the
% boundary of the lid.
% data(i).diastolic.full is the mesh including endo and epi points.

data(i).diastolic.full.xyz = [data(i).diastolic.endo.xyz ; data(i).diastolic.epi.xyz ];
data(i).diastolic.full.tri = [data(i).diastolic.endo.tri ; data(i).diastolic.epi.tri + size(data(i).diastolic.endo.xyz,1) ];
data(i).diastolic.full.B = vtkFeatureEdges( vtkCleanPolyData(data(i).diastolic.full) , 'BoundaryEdgesOn' , [] , 'FeatureEdgesOff' , [] );  %extracting the boundary
data(i).diastolic.full.B.xyz = data(i).diastolic.full.B.xyz( [2 1 3:end], : );  %fixing the connectivity.

data(i).systolic.full.xyz = [data(i).systolic.endo.xyz ; data(i).systolic.epi.xyz ];
data(i).systolic.full.tri = [data(i).systolic.endo.tri ; data(i).systolic.epi.tri + size(data(i).systolic.endo.xyz,1) ];
data(i).systolic.full.B = vtkFeatureEdges( vtkCleanPolyData(data(i).systolic.full) , 'BoundaryEdgesOn' , [] , 'FeatureEdgesOff' , [] );  %extracting the boundary
data(i).systolic.full.B.xyz = data(i).systolic.full.B.xyz( [2 1 3:end], : );  %fixing the connectivity.

%first arg = a mesh, second arg = a list of points
data(i).diastolic.full.B.xyz = data(i).diastolic.full.xyz( vtkClosestPoint( data(i).diastolic.full , data(i).diastolic.full.B.xyz ) , : );
data(i).systolic.full.B.xyz = data(i).systolic.full.xyz( vtkClosestPoint( data(i).systolic.full , data(i).systolic.full.B.xyz ) , : );

% plot3(data(i).diastolic.endo.B.xyz(:,1), data(i).diastolic.endo.B.xyz(:,2), data(i).diastolic.endo.B.xyz(:,3))

%load a struct containing a manually produced myoB.tri
load('C:\Users\jesu2687\Documents\MATLAB\ONBI_short_project_1\myoB_tri.mat')

%join the list of coordinates for endo and epi to be used to cover the myo.
%diastolic
data(i).diastolic.myoB.xyz = [data(i).diastolic.endo.B.xyz ; data(i).diastolic.epi.B.xyz];
data(i).systolic.myoB.xyz = [data(i).systolic.endo.B.xyz ; data(i).systolic.epi.B.xyz];

%append epi surface, endo surface and myo lid
data(i).diastolic.myoB_full.xyz = [ data(i).diastolic.endo.xyz ; data(i).diastolic.myoB.xyz ; data(i).diastolic.epi.xyz ];
data(i).diastolic.myoB_full.tri = [ data(i).diastolic.endo.tri; myoB.tri + size( data(1).diastolic.endo.xyz , 1 ) ; data(i).diastolic.epi.tri + size( data(1).diastolic.endo.xyz , 1 ) + size(data(i).diastolic.myoB.xyz,1) ];
data(i).systolic.myoB_full.xyz = [ data(i).systolic.endo.xyz ; data(i).systolic.myoB.xyz ; data(i).systolic.epi.xyz ];
data(i).systolic.myoB_full.tri = [ data(i).systolic.endo.tri; myoB.tri + size( data(1).systolic.endo.xyz , 1 ) ; data(i).systolic.epi.tri + size( data(1).systolic.endo.xyz , 1 ) + size(data(i).systolic.myoB.xyz,1) ];

%make sure that every triangle points outwards.
data(i).diastolic.myoB_full = FixNormals(data(i).diastolic.myoB_full );
data(i).systolic.myoB_full = FixNormals(data(i).systolic.myoB_full );


%size(transformed_data(i).diastolic.myoB_full)

%calculate volume of myocardium.
[data(i).diastolic.myoVolume, data(i).diastolic.myoCenterOfMass] = MeshVolume( data(i).diastolic.myoB_full );
[data(i).systolic.myoVolume, data(i).systolic.myoCenterOfMass] = MeshVolume( data(i).systolic.myoB_full );
%store volumes
diastolic_myovolumes(i,1) = data(i).diastolic.myoVolume;
systolic_myovolumes(i,1) = data(i).systolic.myoVolume;

% %transformed_data(i).diastolic.myodifference_volume = prod(diff( BBMesh( transformed_data(i).diastolic.myoB_full ) , 1  , 1 ) ) - transformed_data(i).diastolic.myoVolume ;   %%it shoud be positive!!
% transformed_data(i).diastolic.myodifference_volume =  transformed_data(i).diastolic.epi.volume - transformed_data(i).diastolic.myoVolume ;   %%it shoud be positive!!
% transformed_data(i).systolic.myodifference_volume = transformed_data(i).systolic.epi.volume - transformed_data(i).systolic.myoVolume ;   %%it shoud be positive!!

end

hold on
nbins = 100;
histogram(diastolic_myovolumes,nbins)
histogram(systolic_myovolumes,nbins)
legend ' diastolic ' ' systolic'
title 'myocardium volumes'
xlabel 'volume'
ylabel 'frequency'


disp('finished calculating myo volumes')
%% myocardium visualisation
close all
%visualise the endo and epi edge points (labelled), with all the other points from
%endo and epi.
figure
%subplot 221
hold on
plot3(data(i).diastolic.myoB.xyz(:,1),data(i).diastolic.myoB.xyz(:,2), data(i).diastolic.myoB.xyz(:,3))
text(data(i).diastolic.myoB.xyz(:,1) ,  data(i).diastolic.myoB.xyz(:,2) ,  data(i).diastolic.myoB.xyz(:,3) , arrayfun( @(id)sprintf('%d',id) , 1:size(data(i).diastolic.myoB.xyz,1) , 'un',false ) )
plot3(data(i).diastolic.myoB_full.xyz(:,1),data(i).diastolic.myoB_full.xyz(:,2), data(i).diastolic.myoB_full.xyz(:,3),'o')
legend('endo and epi edge points (myo B)', 'all endo and epi points')


%visualise surfaces: endo, epi and myo lid
%cla
figure
%subplot 222
%**************Note: NEED TO CHANGE PROCRUSTES. For myocardium volume calculations (and visualisations) diastolic endo
%and epi should be combined into one long vector. Otherwise, the
%transformation of each will differ.*******************************
patch('vertices',data(i).diastolic.endo.xyz,'faces',data(i).diastolic.endo.tri,'facecolor','none','EdgeColor','red')
patch('vertices',data(i).diastolic.epi.xyz,'faces',data(i).diastolic.epi.tri,'facecolor','none','EdgeColor','blue')
patch('vertices',data(i).diastolic.myoB.xyz,'faces',myoB.tri,'facecolor','green')
legend('endo', 'epi', 'myo lid')

% visualise the myocardium mesh (solid)
%cla
figure
%subplot 223
patch('vertices',data(i).diastolic.myoB_full.xyz,'faces',data(i).diastolic.myoB_full.tri,'facecolor','red')
legend('full myo')

figure
%subplot 223
patch('vertices',data(i).systolic.myoB_full.xyz,'faces',data(i).systolic.myoB_full.tri,'facecolor','red')
legend('full myo')


%visualise center of mass within myocardium mesh
%cla
figure
%subplot 224
patch('vertices',data(i).diastolic.myoB_full.xyz,'faces',data(i).diastolic.myoB_full.tri,'facecolor','g','facealpha',0.1);
hold on;
plot3( data(i).diastolic.myoCenterOfMass(1) ,  data(i).diastolic.myoCenterOfMass(2) ,  data(i).diastolic.myoCenterOfMass(3) , '*r','markers',20 ); hold off
legend('full myo', 'myo center of mass')

%% Shape modeling
disp('start shape modeling')
% PCA
% find mean shape and then use it to find the covariance matrix
[dia_endo_covariance_matrix, dia_endo_mean_shape] = calcCovarianceMatrix(data);

%find eigenvectors of covariance matrix.
%columns of 'dia_endo_eigenvectors' are the eigenvectors.
[dia_endo_eigenvectors, dia_endo_eigenvalues] = eig(dia_endo_covariance_matrix);
%************how do I plot these eigenvectors?...*********
%****should I find covariance of x, y and z seperately?...************

%find the positions of the greatest eigenvalues
[rows, cols] = find((dia_endo_eigenvalues)/max(max(dia_endo_eigenvalues))>=(1))
%select eigenvectors with greatest eigenvalues
dia_endo_eigenvectors(:,1:(cols(1,1)-1)) = 0;

dia_endo_eigenvectors_sum = zeros(size(dia_endo_eigenvectors,2),1);
for i = size(dia_endo_eigenvectors,2)
    dia_endo_eigenvectors_sum = dia_endo_eigenvectors(:,i) + dia_endo_eigenvectors_sum ;
end
dia_endo_new_shape = dia_endo_mean_shape + dia_endo_eigenvectors_sum;

%ICA
%[icaOut] = fastica(dia_endo_covariance_matrix)

dia_endo_new_shape = reshape(dia_endo_new_shape, size(data(1).diastolic.endo.xyz));
disp('finished shape modeling')
%% Shape modelling - visualisation of new shape
hold on
plot3(dia_endo_new_shape(:,1),dia_endo_new_shape(:,2), dia_endo_new_shape(:,3),'.')
dia_endo_mean_shape = reshape(dia_endo_mean_shape, size(data(1).diastolic.endo.xyz));
plot3(dia_endo_mean_shape(:,1),dia_endo_mean_shape(:,2), dia_endo_mean_shape(:,3),'+')

%%

% phi = dia_endo_eigenvectors';
% phi = (sortrows(phi))';

%reverse the eigenvalue matrix
dia_endo_eigenvalues = dia_endo_eigenvalues(end:-1:1);
% reverse the columns
dia_endo_eigenvectors = dia_endo_eigenvectors(:,end:-1:1); 
dia_endo_eigenvectors=dia_endo_eigenvectors';

% model parameters 'b'

%phi = dia_endo_eigenvectors(:);

% select some of the largest eigenvalues
%eigenvalue_indices = find(phi>0.6, 20);
%phi(eigenvalue_indices);

% sortedPhi = sort(phi, 'descend');
% phi = reshape((sortedPhi.'), size(dia_endo_covariance_matrix)); %each column shows an eigenvetor

% %coeff = pca(covariance_matrix)

% calculate new shapes, using different parameters b
%!!!!!!!!!!!WHAT VALUES SHOULD 'b' BE?
% new_shape = mean_shape + (phi)*(b);
%%
for i = 1:size(dia_endo_eigenvectors,2)
b(1,i) = (dia_endo_eigenvectors(:,i).')*(data(1).diastolic.endo.xyz(:) - dia_endo_mean_shape); %transpose phi
end
%%
% plot in pca space... not sure if this is correct...
for i = 1:400
pc1 = eigenvectors * data(i).diastolic.endo.xyz;
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
%!!!!! NEED TO CORRECT THIS...ALREADY DONE BY 'calcVolume'?
endo_centroid = mean(endo);
epi_centroid = mean(epi);
