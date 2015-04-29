
%function[ myoVolumes ] = calcMyoVolumes(

setenv('path',[getenv('path'),';','C:\Users\jesu2687\Documents\MATLAB\ErnestoCode\Tools\MESHES\vtk_libs']);
addpath C:\Users\jesu2687\Documents\MATLAB\ErnestoCode\Tools\MESHES\
addpath C:\Users\jesu2687\Documents\MATLAB\ErnestoCode\Tools\
 
cd C:\Users\jesu2687\Documents\MATLAB\ErnestoCode\


 % vtkCleanPolyData(EPI_ED) fix the possible replicated nodes and spurious
 % edges
test.diastolic.endo.B = vtkFeatureEdges( vtkCleanPolyData(transformed_data(1).diastolic.endo) , 'BoundaryEdgesOn' , [] , 'FeatureEdgesOff' , [] );  %extracting the boundary
test.diastolic.endo.B.xyz = test.diastolic.endo.B.xyz( [2 1 3:end], : );  %fixing the connectivity.

test.diastolic.epi.B = vtkFeatureEdges( vtkCleanPolyData(transformed_data(1).diastolic.epi) , 'BoundaryEdgesOn' , [] , 'FeatureEdgesOff' , [] );  %extracting the boundary
test.diastolic.epi.B.xyz = test.diastolic.epi.B.xyz( [2 1 3:end], : );  %fixing the connectivity.

%load a manually produced myoB.tri
load('C:\Users\jesu2687\Documents\MATLAB\ONBI_short_project_1\myoB.mat')

%join the list of coordinates for endo and epi to be used to cover the myo
%  myoB.endo.xyz = test.diastolic.endo.B.xyz;
%  myoB.epi.xyz = test.diastolic.epi.B.xyz;
%  myoB.xyz = [myoB.endo.xyz ; myoB.epi.xyz];
 
%myoB.xyz stored in myoB.mat
plot3(myoB.xyz(:,1),myoB.xyz(:,2), myoB.xyz(:,3))
% visualise labelled points on the edge of the myocardium lid
plot3(myoB.xyz(:,1) ,  myoB.xyz(:,2) ,  myoB.xyz(:,3) , '.r' )
text(myoB.xyz(:,1) ,  myoB.xyz(:,2) ,  myoB.xyz(:,3) , arrayfun( @(id)sprintf('%d',id) , 1:size(myoB.xyz,1) , 'un',false ) )


% %visualise surfaces: endo, epi and myo lid
% cla
% patch('vertices',data(1).diastolic.endo.xyz,'faces',data(1).diastolic.endo.tri,'facecolor','none','EdgeColor','red')
% patch('vertices',data(1).diastolic.epi.xyz,'faces',data(1).diastolic.epi.tri,'facecolor','none','EdgeColor','blue')
% patch('vertices',myoB.xyz,'faces',myoB.tri,'facecolor','g')

%append epi surface, endo surface and myo lid
% Last_prevID  = size( EPI_ED.xyz , 1 );
% EPI_ED.xyz = [ EPI_ED.xyz ; B.xyz ];
% EPI_ED.tri = [ EPI_ED.tri ; B.tri + Last_prevID ];

Last_prevID1  = size( data(1).diastolic.endo.xyz , 1 );
Last_prevID2  = size( data(1).diastolic.epi.xyz , 1 );
Last_prevID3 = size(myoB.xyz,1);
myoB_full.xyz = [ data(1).diastolic.endo.xyz ; myoB.xyz ; data(1).diastolic.epi.xyz ];
myoB_full.tri = [ data(1).diastolic.endo.tri; myoB.tri + Last_prevID1 ; data(1).diastolic.epi.tri + Last_prevID1 + Last_prevID3 ];

%make sure that every triangle points outwards
myoB_full = FixNormals( myoB_full );
size(myoB_full)
% visualise myocardium mesh

figure
patch('vertices',myoB_full.xyz,'faces',myoB_full.tri,'facecolor','red')

% hold on
% plot3(myoB.xyz(:,1),myoB.xyz(:,2), myoB.xyz(:,3),'+')
% plot3(myoB_full.xyz(:,1),myoB_full.xyz(:,2), myoB_full.xyz(:,3),'o')

%calculate volume of myocardium
[Volume,CenterOfMass] = MeshVolume( myoB_full )

difference_volume = prod(diff( BBMesh( myoB_full ) , 1  , 1 ) ) - Volume   %%it shoud be positive!!

%visualise center of mass within myocardium mesh
cla
patch('vertices',myoB_full.xyz,'faces',myoB_full.tri,'facecolor','g','facealpha',0.1);
hold on;
plot3( CenterOfMass(1) , CenterOfMass(2) , CenterOfMass(3) , '*r','markers',20 ); hold off







