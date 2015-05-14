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
% find nearest points on endo and epi boundaries and assign them as the
% boundary of the lid.
% data(i).diastolic.full is the mesh including endo and epi points.
data(i).diastolic.full.xyz = [data(i).diastolic.endo.xyz ; data(i).diastolic.epi.xyz ];
data(i).diastolic.full.tri = [data(i).diastolic.endo.tri ; data(i).diastolic.epi.tri + size(data(i).diastolic.endo.xyz,1) ];
data(i).diastolic.full.B = vtkFeatureEdges( vtkCleanPolyData(data(i).diastolic.full) , 'BoundaryEdgesOn' , [] , 'FeatureEdgesOff' , [] );  %extracting the boundary
data(i).diastolic.full.B.xyz = data(i).diastolic.full.B.xyz( [2 1 3:end], : );  %fixing the connectivity.

%first arg = a mesh, second arg = a list of points
data(i).diastolic.full.B.xyz = data(i).diastolic.full.xyz( vtkClosestPoint( data(i).diastolic.full , data(i).diastolic.full.B.xyz ) , : );

% plot3(data(i).diastolic.endo.B.xyz(:,1), data(i).diastolic.endo.B.xyz(:,2), data(i).diastolic.endo.B.xyz(:,3))

%load a struct containing a manually produced myoB.tri
load('C:\Users\jesu2687\Documents\MATLAB\ONBI_short_project_1\myoB_tri.mat')

%join the list of coordinates for endo and epi to be used to cover the myo.
%diastolic
data(i).diastolic.myoB.xyz = [data(i).diastolic.endo.B.xyz ; data(i).diastolic.epi.B.xyz];

%append epi surface, endo surface and myo lid
%diastolic

data(i).diastolic.myoB_full.xyz = [ data(i).diastolic.endo.xyz ; data(i).diastolic.myoB.xyz ; data(i).diastolic.epi.xyz ];
data(i).diastolic.myoB_full.tri = [ data(i).diastolic.endo.tri; myoB.tri + size( data(1).diastolic.endo.xyz , 1 ) ; data(i).diastolic.epi.tri + size( data(1).diastolic.endo.xyz , 1 ) + size(data(i).diastolic.myoB.xyz,1) ];

%make sure that every triangle points outwards.
%diastolic
data(i).diastolic.myoB_full = FixNormals(data(i).diastolic.myoB_full );


%size(transformed_data(i).diastolic.myoB_full)

%calculate volume of myocardium.
%diastolic
[data(i).diastolic.myoVolume, data(i).diastolic.myoCenterOfMass] = MeshVolume( data(i).diastolic.myoB_full );

% %transformed_data(i).diastolic.myodifference_volume = prod(diff( BBMesh( transformed_data(i).diastolic.myoB_full ) , 1  , 1 ) ) - transformed_data(i).diastolic.myoVolume ;   %%it shoud be positive!!
% transformed_data(i).diastolic.myodifference_volume =  transformed_data(i).diastolic.epi.volume - transformed_data(i).diastolic.myoVolume ;   %%it shoud be positive!!
% transformed_data(i).systolic.myodifference_volume = transformed_data(i).systolic.epi.volume - transformed_data(i).systolic.myoVolume ;   %%it shoud be positive!!

end

disp('finished calculating myo volumes')