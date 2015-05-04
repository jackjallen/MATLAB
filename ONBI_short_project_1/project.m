%% ONBI Short project 1
% Jack Allen
% Supervisor: Vicente Grau
%
clear all
close all
clc

<<<<<<< HEAD
setenv('path',[getenv('path'),';','C:\Users\jesu2687\Documents\MATLAB\ErnestoCode\Tools\MESHES\vtk_libs']);
addpath('C:\Users\jesu2687\Documents\MATLAB\ONBI_short_project_1')
addpath C:\Users\jesu2687\Documents\MATLAB\ErnestoCode\Tools\MESHES\
addpath C:\Users\jesu2687\Documents\MATLAB\ErnestoCode\Tools
addpath C:\Users\jesu2687\Documents\MATLAB\ErnestoCode\
load('C:\Users\jesu2687\Documents\MATLAB\ONBI_short_project_1\project1_data.mat');
load('C:\Users\jesu2687\Documents\MATLAB\ONBI_short_project_1\labels.mat');

%store indices that allow us to indentify which study (DETERMINE or MESA)
%each case belongs to.
% (DETERMINE = mycardium has been infarcted)
% (MESA = asymptomatic)
=======
% setenv('path',[getenv('path'),';','C:\Users\jesu2687\Documents\MATLAB\ErnestoCode\Tools\MESHaES\vtk_libs']);
% addpath('C:\Users\jesu2687\Documents\MATLAB\ONBI_short_project_1')
% addpath C:\Users\jesu2687\Documents\MATLAB\ErnestoCode\Tools\MESHES\
% addpath C:\Users\jesu2687\Documents\MATLAB\ErnestoCode\Tools
% addpath C:\Users\jesu2687\Documents\MATLAB\ErnestoCode\
% 
setenv('path',[getenv('path'),';','D:\ErnestoCode\Tools\MESHES\vtk_libs']);
addpath('C:\Users\jack laptop\Documents\MATLAB\ONBI_DTC\short_project_1')
addpath ('C:\Users\jack laptop\Documents\MATLAB\ONBI_DTC\short_project_1\ErnestoCode\Tools\MESHES\');
addpath ('C:\Users\jack laptop\Documents\MATLAB\ONBI_DTC\short_project_1\ErnestoCode\Tools');
addpath ('C:\Users\jack laptop\Documents\MATLAB\ONBI_DTC\short_project_1\ErnestoCode');
% read data (adapted from VG code)
%load('C:\Users\jesu2687\Documents\MATLAB\ONBI_short_project_1\project1_data.mat');
load('C:\Users\jack laptop\Documents\MATLAB\ONBI_DTC\short_project_1\project1_data.mat');
load('C:\Users\jack laptop\Documents\MATLAB\ONBI_DTC\short_project_1\labels.mat');

>>>>>>> 5ebf88f09de8911523cf4d9c2ee06a809999e990
DETERMINE_indices = Labels(2:101,3);
MESA_indices = Labels(102:201,3);

%% Procrustes analysis
% This is done by concatenating the endo and epi shape vectors to produce
% longer vectors (myo.xyz). Procrustes analysis is then performed on the
% two myo.xyz vectors. Finally, the transformed myo.xyz vectors are split
% to extract transformed endo and epi shape vectors.
disp('starting procrustes analysis')

%store dimensions of the shapes
[shape_nRows , shape_nCols] = size(data(1).diastolic.endo.xyz);

%use the initial data clouds as references (individual)
dia_endo_reference = data(1).diastolic.endo.xyz;
dia_epi_reference = data(1).diastolic.epi.xyz;
sys_endo_reference = data(1).systolic.endo.xyz;
sys_epi_reference = data(1).systolic.epi.xyz;
%make a reference from the endo and epi surfaces (concatenate)
dia_myo_reference = [dia_endo_reference(:) ; dia_epi_reference(:)];
sys_myo_reference = [sys_endo_reference(:)  ; sys_epi_reference(:)];

%initialise the sums of the shapes as zero
% dia_endo_sum = zeros(size(dia_endo_reference));
% dia_epi_sum = zeros(size(dia_epi_reference));
% sys_endo_sum = zeros(size(sys_endo_reference));
% sys_epi_sum = zeros(size(sys_epi_reference));
% (whole shape)
dia_myo_sum = zeros(size(dia_myo_reference));
sys_myo_sum = zeros(size(sys_myo_reference));

% Think of endo and epi as one shape (concatenate).
for i = 1:400
data(i).diastolic.myo.xyz = [data(i).diastolic.endo.xyz ; data(i).diastolic.epi.xyz];
data(i).systolic.myo.xyz = [data(i).systolic.endo.xyz ; data(i).systolic.epi.xyz];
end

%p = number of times procrustes is performed.
for p = 1:2
% iterate so that procrustes is performed on each case (400 patients).
for i = 1:400
% % transform endo and epi (individually)
% [data(i).diastolic.endo.procrustes_d, data(i).diastolic.endo.xyz] = procrustes(dia_endo_reference, data(i).diastolic.endo.xyz(:));
% [data(i).diastolic.epi.procrustes_d, data(i).diastolic.epi.xyz] = procrustes(dia_epi_reference, data(i).diastolic.epi.xyz(:));
% [data(i).systolic.endo.procrustes_d, data(i).systolic.endo.xyz] = procrustes(sys_endo_reference, data(i).systolic.endo.xyz(:));
% [data(i).systolic.epi.procrustes_d, data(i).systolic.epi.xyz] = procrustes(sys_epi_reference, data(i).systolic.epi.xyz(:));

% transform endo and epi as one shape vector
[data(i).diastolic.myo.procrustes_d, data(i).diastolic.myo.xyz] = procrustes(dia_myo_reference, data(i).diastolic.myo.xyz(:));
[data(i).systolic.myo.procrustes_d, data(i).systolic.myo.xyz] = procrustes(sys_myo_reference, data(i).systolic.myo.xyz(:));

%sums for finding new means later on
% dia_endo_sum = data(i).diastolic.endo.xyz + dia_endo_sum;
% dia_epi_sum = data(i).diastolic.epi.xyz + dia_epi_sum;
% sys_endo_sum = data(i).systolic.endo.xyz + sys_endo_sum;
% sys_epi_sum = data(i).systolic.epi.xyz + sys_epi_sum;
dia_myo_sum = data(i).diastolic.myo.xyz + dia_myo_sum;
sys_myo_sum = data(i).systolic.myo.xyz + sys_myo_sum;

end

%calculate means
% dia_endo_mean = dia_endo_sum/400;
% dia_epi_mean = dia_epi_sum/400;
% sys_endo_mean = sys_endo_sum/400;
% sys_epi_mean = sys_epi_sum/400;
dia_myo_mean = dia_myo_sum/400;
sys_myo_mean = sys_myo_sum/400;

%set means as the references for next round of procrustes
% dia_endo_reference = dia_endo_mean;
% dia_epi_reference = dia_epi_mean;
% sys_endo_reference = sys_endo_mean;
% sys_epi_reference = sys_epi_mean;
dia_myo_reference = dia_myo_mean; 
sys_myo_reference = sys_myo_mean;
end

% % reshape from vector to matrix.
% % individual endo and epi shapes.
% dia_endo_mean = reshape(dia_endo_mean,size(data(1).diastolic.endo.xyz));
% dia_epi_mean = reshape(dia_epi_mean,size(data(1).diastolic.endo.xyz));
% sys_endo_mean = reshape(sys_endo_mean,size(data(1).diastolic.endo.xyz));
% sys_epi_mean = reshape(sys_epi_mean,size(data(1).diastolic.endo.xyz));

% % concatenated endo and epi shapes.
% % split the long vectors into vectors that represent the endo and epi
% % shapes then reshape the resulting vectors to matrix form.
dia_myo_reference = reshape(dia_myo_reference, [2*shape_nRows, shape_nCols]);
sys_myo_reference = reshape(sys_myo_reference, [2*shape_nRows, shape_nCols]);
for i = 1:400
% reshape myo shapes from vector to matrix
data(i).diastolic.myo.xyz = reshape(data(i).diastolic.myo.xyz,[2*shape_nRows, shape_nCols]);
data(i).systolic.myo.xyz = reshape(data(i).systolic.myo.xyz,[2*shape_nRows, shape_nCols]);

% extract endo and epi shapes from full shape matrices (data(i).diastolic.myo.xyz)
data(i).diastolic.endo.xyz = data(i).diastolic.myo.xyz(1:shape_nRows, :);
data(i).diastolic.epi.xyz = data(i).diastolic.myo.xyz(shape_nRows+1:2*shape_nRows, :);
data(i).systolic.endo.xyz = data(i).systolic.myo.xyz(1:shape_nRows, :);
data(i).systolic.epi.xyz = data(i).systolic.myo.xyz(shape_nRows+1:2*shape_nRows, :);
end
disp('finished procrustes analysis')
%% Visualising to check procrustes output
disp('visualise procrustes output')

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
hold on
plot3(data(i).systolic.myo.xyz(:,1), data(i).systolic.myo.xyz(:,2), data(i).systolic.myo.xyz(:,3),'o'); 
title 'patient i, whole LV, systolic'

%% make lid for myocardium
% disp('make lid for myocardium')
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

disp('finish calculating endo and epi volumes')
<<<<<<< HEAD
%% plot volume histograms
=======
%% Myocardium volumes
disp('started calculating myocardium volumes')
% for i = 1:400 %all patients
load('myoB.mat')


>>>>>>> 5ebf88f09de8911523cf4d9c2ee06a809999e990




%% Myocardium volumes
disp('calculating myocardium volumes')

%load struct containing the triangle file for the myo 'donut' shaped lid
load('C:\Users\jesu2687\Documents\MATLAB\ONBI_short_project_1\myoB_tri.mat') 

% for i = 1:400 %all patients
for i = 1:400;
% Find endo and epi boundary points (B.xyz)
% vtkCleanPolyData(EPI_ED) fix the possible replicated nodes and spurious
% edges.
% data(i).diastolic.endo.B = vtkFeatureEdges( vtkCleanPolyData(data(i).diastolic.endo) , 'BoundaryEdgesOn' , [] , 'FeatureEdgesOff' , [] );  %extracting the boundary
% data(i).diastolic.endo.B.xyz = data(i).diastolic.endo.B.xyz( [2 1 3:end], : );  %fixing the connectivity.
% data(i).diastolic.epi.B = vtkFeatureEdges( vtkCleanPolyData(data(i).diastolic.epi) , 'BoundaryEdgesOn' , [] , 'FeatureEdgesOff' , [] );  %extracting the boundary
% data(i).diastolic.epi.B.xyz = data(i).diastolic.epi.B.xyz( [2 1 3:end], : );  %fixing the connectivity.
% data(i).systolic.endo.B = vtkFeatureEdges( vtkCleanPolyData(data(i).systolic.endo) , 'BoundaryEdgesOn' , [] , 'FeatureEdgesOff' , [] );  %extracting the boundary
% data(i).systolic.endo.B.xyz = data(i).systolic.endo.B.xyz( [2 1 3:end], : );  %fixing the connectivity.
% data(i).systolic.epi.B = vtkFeatureEdges( vtkCleanPolyData(data(i).systolic.epi) , 'BoundaryEdgesOn' , [] , 'FeatureEdgesOff' , [] );  %extracting the boundary
% data(i).systolic.epi.B.xyz = data(i).systolic.epi.B.xyz( [2 1 3:end], : );  %fixing the connectivity.

%Find full shape boundary points (B.xyz)
% full = full shape without myo lid
data(i).diastolic.full.xyz = [data(i).diastolic.endo.xyz ; data(i).diastolic.epi.xyz ];
data(i).diastolic.full.tri = [data(i).diastolic.endo.tri ; data(i).diastolic.epi.tri + size(data(i).diastolic.endo.xyz,1) ];
data(i).diastolic.full.B = vtkFeatureEdges( vtkCleanPolyData(data(i).diastolic.full) , 'BoundaryEdgesOn' , [] , 'FeatureEdgesOff' , [] );  %extracting the boundary
data(i).diastolic.full.B.xyz = data(i).diastolic.full.B.xyz( [2 1 3:end], : );  %fixing the connectivity.
data(i).systolic.full.xyz = [data(i).systolic.endo.xyz ; data(i).systolic.epi.xyz ];
data(i).systolic.full.tri = [data(i).systolic.endo.tri ; data(i).systolic.epi.tri + size(data(i).systolic.endo.xyz,1) ];
data(i).systolic.full.B = vtkFeatureEdges( vtkCleanPolyData(data(i).systolic.full) , 'BoundaryEdgesOn' , [] , 'FeatureEdgesOff' , [] );  %extracting the boundary
data(i).systolic.full.B.xyz = data(i).systolic.full.B.xyz( [2 1 3:end], : );  %fixing the connectivity.
% find nearest points on endo and epi boundaries (identically positioned, but not connected to main shape) and assign them as the
% boundary of the lid.
% first arg = a mesh, second arg = a list of point coordinates.
data(i).diastolic.full.B.xyz = data(i).diastolic.full.xyz( vtkClosestPoint( data(i).diastolic.full , data(i).diastolic.full.B.xyz ) , : );
data(i).systolic.full.B.xyz = data(i).systolic.full.xyz( vtkClosestPoint( data(i).systolic.full , data(i).systolic.full.B.xyz ) , : );

<<<<<<< HEAD
%join the list of coordinates for endo and epi to be used to make myo lid
data(i).diastolic.myo.B.xyz = [data(i).diastolic.endo.B.xyz ; data(i).diastolic.epi.B.xyz]; 
data(i).systolic.myo.B.xyz = [data(i).systolic.endo.B.xyz ; data(i).systolic.epi.B.xyz]; 
=======
% plot3(data(i).diastolic.endo.B.xyz(:,1), data(i).diastolic.endo.B.xyz(:,2), data(i).diastolic.endo.B.xyz(:,3))

%load a struct containing a manually produced myoB.tri
% load('C:\Users\jesu2687\Documents\MATLAB\ONBI_short_project_1\myoB_tri.mat')

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
>>>>>>> 5ebf88f09de8911523cf4d9c2ee06a809999e990

%Make fully closed myocardium volume by appending epi surface, endo
%surface and myo lid.
data(i).diastolic.myo.xyz = [ data(i).diastolic.endo.xyz ; data(i).diastolic.myo.B.xyz ; data(i).diastolic.epi.xyz ];
data(i).diastolic.myo.tri = [ data(i).diastolic.endo.tri; myoB.tri + size( data(1).diastolic.endo.xyz , 1 ) ; data(i).diastolic.epi.tri + size( data(1).diastolic.endo.xyz , 1 ) + size(data(i).diastolic.myo.B.xyz,1) ];
data(i).systolic.myo.xyz = [ data(i).systolic.endo.xyz ; data(i).systolic.myo.B.xyz ; data(i).systolic.epi.xyz ];
data(i).systolic.myo.tri = [ data(i).systolic.endo.tri; myoB.tri + size( data(1).systolic.endo.xyz , 1 ) ; data(i).systolic.epi.tri + size( data(1).systolic.endo.xyz , 1 ) + size(data(i).systolic.myo.B.xyz,1) ];

%make sure that the normal to each triangle points outwards.
data(i).diastolic.myo = FixNormals(data(i).diastolic.myo );
data(i).systolic.myo = FixNormals(data(i).systolic.myo );

%calculate volume of myocardium.
<<<<<<< HEAD
[data(i).diastolic.myoVolume, data(i).diastolic.myoCenterOfMass] = MeshVolume( data(i).diastolic.myo );
[data(i).systolic.myoVolume, data(i).systolic.myoCenterOfMass] = MeshVolume( data(i).systolic.myo );
=======
[data(i).diastolic.myoVolume, data(i).diastolic.myoCenterOfMass] = MeshVolume( data(i).diastolic.myoB_full );
[data(i).systolic.myoVolume, data(i).systolic.myoCenterOfMass] = MeshVolume( data(i).systolic.myoB_full );
>>>>>>> 5ebf88f09de8911523cf4d9c2ee06a809999e990

% %Compare myo volume with epi volume.
% %transformed_data(i).diastolic.myodifference_volume = prod(diff( BBMesh( transformed_data(i).diastolic.myoB_full ) , 1  , 1 ) ) - transformed_data(i).diastolic.myoVolume ;   %%it shoud be positive!!
% transformed_data(i).diastolic.myodifference_volume =  transformed_data(i).diastolic.epi.volume - transformed_data(i).diastolic.myoVolume ;   %%it shoud be positive!!
% transformed_data(i).systolic.myodifference_volume = transformed_data(i).systolic.epi.volume - transformed_data(i).systolic.myoVolume ;   %%it shoud be positive!!

end

<<<<<<< HEAD
disp('finished calculating myocardium volumes')
=======
>>>>>>> 5ebf88f09de8911523cf4d9c2ee06a809999e990
%% store volumes
DETERMINE_indices = sort(DETERMINE_indices);
for i = DETERMINE_indices(1,1):DETERMINE_indices(100,1)
DETERMINE_diastolic_myovolumes(i,1) = data(i).diastolic.myoVolume;
DETERMINE_systolic_myovolumes(i,1) = data(i).systolic.myoVolume;
end
MESA_indices = sort(MESA_indices);
for i = MESA_indices(1,1):MESA_indices(100,1)
MESA_diastolic_myovolumes(i,1) = data(i).diastolic.myoVolume;
MESA_systolic_myovolumes(i,1) = data(i).systolic.myoVolume;
end

figure
<<<<<<< HEAD
hold on
nbins = 100;
histogram(DETERMINE_systolic_myovolumes,nbins)
histogram(MESA_systolic_myovolumes,nbins)
legend 'DETERMINE' 'MESA'
title 'systolic myocardium volumes'
=======
hold on
nbins = 100;
histogram(DETERMINE_systolic_myovolumes,nbins)
histogram(MESA_systolic_myovolumes,nbins)
legend 'DETERMINE' 'MESA'
title 'systolic myocardium volumes'
xlabel 'volume'
ylabel 'frequency'

figure
hold on
nbins = 100;
histogram(DETERMINE_diastolic_myovolumes,nbins)
histogram(MESA_diastolic_myovolumes,nbins)
legend ' DETERMINE ' ' MESA '
title 'diastolic myocardium volumes'
>>>>>>> 5ebf88f09de8911523cf4d9c2ee06a809999e990
xlabel 'volume'
ylabel 'frequency'

figure
hold on
nbins = 100;
histogram(DETERMINE_diastolic_myovolumes,nbins)
histogram(MESA_diastolic_myovolumes,nbins)
legend ' DETERMINE ' ' MESA '
title 'diastolic myocardium volumes'
xlabel 'volume'
ylabel 'frequency'

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
