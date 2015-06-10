
setenv('path',[getenv('path'),';','F:\ErnestoCode\Tools\MESHES\vtk_libs']);
addpath F:\ErnestoCode\Tools\MESHES\
addpath F:\ErnestoCode\Tools\

cd F:\ErnestoCode\

 % vtkCleanPolyData(EPI_ED) fix the possible replicated nodes and spurious
 % edges
test.diastolic.endo.B = vtkFeatureEdges( vtkCleanPolyData(transformed_data(1).diastolic.endo) , 'BoundaryEdgesOn' , [] , 'FeatureEdgesOff' , [] );  %extracting the boundary
test.diastolic.endo.B.xyz = test.diastolic.endo.B.xyz( [2 1 3:end], : );  %fixing the connectivity.

test.diastolic.epi.B = vtkFeatureEdges( vtkCleanPolyData(transformed_data(1).diastolic.epi) , 'BoundaryEdgesOn' , [] , 'FeatureEdgesOff' , [] );  %extracting the boundary
test.diastolic.epi.B.xyz = test.diastolic.epi.B.xyz( [2 1 3:end], : );  %fixing the connectivity.

%load a manually produced myoB.tri
load('C:\Users\jesu2687\Documents\MATLAB\ONBI_short_project_1\myoB.mat')

%join the list of coordinates for endo and epi to be used to cover the myo
 myoB.endo.xyz = test.diastolic.endo.B.xyz;
 myoB.epi.xyz = test.diastolic.epi.B.xyz;
 myoB.xyz = [myoB.epi.xyz ; myoB.endo.xyz];
 
 plot3(myoB.xyz(:,1),myoB.xyz(:,2), myoB.xyz(:,3))
 
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
myoB_full.xyz = [ data(1).diastolic.endo.xyz ; myoB.xyz; data(1).diastolic.epi.xyz ];
myoB_full.tri = [ data(1).diastolic.endo.tri; myoB.tri + Last_prevID1 ; data(1).diastolic.epi.tri ];

%make sure that every triangle points outwards
myoB_full = FixNormals( myoB_full );

cla
figure
patch('vertices',myoB_full.xyz,'faces',myoB_full.tri,'facecolor','r')

plot3(myoB.xyz(:,1),myoB.xyz(:,2), myoB.xyz(:,3),'o')
plot3(myoB_full.xyz(:,1),myoB_full.xyz(:,2), myoB_full.xyz(:,3),'o')

[Volume,CenterOfMass] = MeshVolume( myoB_full )

difference_volume = prod(diff( BBMesh( EPI_ED ) , 1  , 1 ) ) - Volume   %%it shoud be positive!!

cla
patch('vertices',EPI_ED.xyz,'faces',EPI_ED.tri,'facecolor','b','facealpha',0.2); hold on; plot3( CenterOfMass(1) , CenterOfMass(2) , CenterOfMass(3) , '*r' ); hold off

%%

plot3( B.xyz(:,1) ,  B.xyz(:,2) ,  B.xyz(:,3) , '.r' )
text(  B.xyz(:,1) ,  B.xyz(:,2) ,  B.xyz(:,3) , arrayfun( @(id)sprintf('%d',id) , 1:size(B.xyz,1) , 'un',false ) )







