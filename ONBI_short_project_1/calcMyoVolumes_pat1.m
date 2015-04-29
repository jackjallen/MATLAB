
 setenv('path',[getenv('path'),';','C:\Users\jesu2687\Documents\MATLAB\ErnestoCode\Tools\MESHES\vtk_libs']);
 addpath C:\Users\jesu2687\Documents\MATLAB\ErnestoCode\Tools\MESHES\
 addpath C:\Users\jesu2687\Documents\MATLAB\ErnestoCode\Tools\
 
cd C:\Users\jesu2687\Documents\MATLAB\ErnestoCode\

for i = 1:400

 % vtkCleanPolyData(EPI_ED) fix the possible replicated nodes and spurious
 % edges
 
 %reshape from a vector to a matrix
transformed_data(i).diastolic.endo.xyz = reshape(transformed_data(i).diastolic.endo.xyz, size(data(i).diastolic.endo.xyz));
transformed_data(i).diastolic.epi.xyz = reshape(transformed_data(i).diastolic.epi.xyz, size(data(i).diastolic.epi.xyz));

transformed_data(i).diastolic.endo.B = vtkFeatureEdges( vtkCleanPolyData(transformed_data(i).diastolic.endo) , 'BoundaryEdgesOn' , [] , 'FeatureEdgesOff' , [] );  %extracting the boundary
transformed_data(i).diastolic.endo.B.xyz = transformed_data(i).diastolic.endo.B.xyz( [2 1 3:end], : );  %fixing the connectivity.

transformed_data(i).diastolic.epi.B = vtkFeatureEdges( vtkCleanPolyData(transformed_data(i).diastolic.epi) , 'BoundaryEdgesOn' , [] , 'FeatureEdgesOff' , [] );  %extracting the boundary
transformed_data(i).diastolic.epi.B.xyz = transformed_data(i).diastolic.epi.B.xyz( [2 1 3:end], : );  %fixing the connectivity.

%load a manually produced myoB.tri
load('C:\Users\jesu2687\Documents\MATLAB\ONBI_short_project_1\myoB.mat')

%join the list of coordinates for endo and epi to be used to cover the myo
 transformed_data(i).diastolic.myoB.xyz = [transformed_data(i).diastolic.endo.B.xyz ; transformed_data(i).diastolic.epi.B.xyz];
 
%append epi surface, endo surface and myo lid
% Last_prevID  = size( EPI_ED.xyz , 1 );
% EPI_ED.xyz = [ EPI_ED.xyz ; B.xyz ];
% EPI_ED.tri = [ EPI_ED.tri ; B.tri + Last_prevID ];

Last_prevID1  = size( data(1).diastolic.endo.xyz , 1 );
Last_prevID2  = size( data(1).diastolic.epi.xyz , 1 );
Last_prevID3 = size(transformed_data(i).diastolic.myoB.xyz,1);
transformed_data(i).diastolic.myoB_full.xyz = [ transformed_data(i).diastolic.endo.xyz ; transformed_data(i).diastolic.myoB.xyz ; transformed_data(i).diastolic.epi.xyz ];
transformed_data(i).diastolic.myoB_full.tri = [ transformed_data(i).diastolic.endo.tri; myoB.tri + Last_prevID1 ; transformed_data(i).diastolic.epi.tri + Last_prevID1 + Last_prevID3 ];

%make sure that every triangle points outwards
transformed_data(i).diastolic.myoB_full = FixNormals(transformed_data(i).diastolic.myoB_full);
% size(transformed_data(i).diastolic.myoB_full)

% hold on
% plot3(myoB.xyz(:,1),myoB.xyz(:,2), myoB.xyz(:,3),'+')
% plot3(myoB_full.xyz(:,1),myoB_full.xyz(:,2), myoB_full.xyz(:,3),'o')

%calculate volume of myocardium
[transformed_data(i).diastolic.myoVolume, transformed_data(i).diastolic.myoCenterOfMass] = MeshVolume( transformed_data(i).diastolic.myoB_full );

%transformed_data(i).diastolic.myodifference_volume = prod(diff( BBMesh( transformed_data(i).diastolic.myoB_full ) , 1  , 1 ) ) - transformed_data(i).diastolic.myoVolume ;   %%it shoud be positive!!
transformed_data(i).diastolic.myodifference_volume = prod(diff( BBMesh( transformed_data(i).diastolic.myoB_full ) , 1  , 1 ) ) - transformed_data(i).diastolic.myoVolume ;   %%it shoud be positive!!
end


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




