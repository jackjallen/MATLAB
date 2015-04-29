
 setenv('path',[getenv('path'),';','C:\Users\jesu2687\Documents\MATLAB\ErnestoCode\Tools\MESHES\vtk_libs']);
 addpath C:\Users\jesu2687\Documents\MATLAB\ErnestoCode\Tools\MESHES\
 addpath C:\Users\jesu2687\Documents\MATLAB\ErnestoCode\Tools\
 
cd C:\Users\jesu2687\Documents\MATLAB\ErnestoCode\

 % vtkCleanPolyData(EPI_ED) fix the possible replicated nodes and spurious
 % edges
for i = 1:400
     
transformed_data(i).diastolic.endo.B = vtkFeatureEdges( vtkCleanPolyData(transformed_data(i).diastolic.endo) , 'BoundaryEdgesOn' , [] , 'FeatureEdgesOff' , [] );  %extracting the boundary
transformed_data(i).diastolic.endo.B.xyz = test.diastolic.endo.B.xyz( [2 1 3:end], : );  %fixing the connectivity.

transformed_data(i).diastolic.epi.B = vtkFeatureEdges( vtkCleanPolyData(transformed_data(i).diastolic.epi) , 'BoundaryEdgesOn' , [] , 'FeatureEdgesOff' , [] );  %extracting the boundary
transformed_data(i).diastolic.epi.B.xyz = test.diastolic.epi.B.xyz( [2 1 3:end], : );  %fixing the connectivity.

%load a manually produced myoB.tri
load('C:\Users\jesu2687\Documents\MATLAB\ONBI_short_project_1\myoB.mat')

%join the list of coordinates for endo and epi to be used to cover the myo
  myoB(i).dia.endo.xyz = transformed_data(i).diastolic.endo.B.xyz;
  myoB(i).dia.epi.xyz = transformed_data(i).diastolic.epi.B.xyz;
  myoB(i).sys.endo.xyz = transformed_data(i).systolic.endo.B.xyz;
  myoB(i).sys.epi.xyz = transformed_data(i).systolic.epi.B.xyz;
  
  myoB(i).dia.xyz = [myoB(i).dia.endo.xyz ; myoB(i).dia.epi.xyz];
  myoB(i).sys.xyz = [myoB(i).sys.endo.xyz ; myoB(i).sys.epi.xyz];

 

plot3(myoB(1).dia.xyz(:,1),myoB(1).dia.xyz(:,2), myoB(1).dia.xyz(:,3))
% visualise labelled points on the edge of the myocardium lid
plot3(myoB(1).dia.xyz(:,1) ,  myoB(1).dia.xyz(:,2) ,  myoB(1).dia.xyz(:,3) , '.r' )
text(myoB(1).dia.xyz(:,1) ,  myoB(1).dia.xyz(:,2) ,  myoB(1).dia.xyz(:,3) , arrayfun( @(id)sprintf('%d',id) , 1:size(myoB(1).dia.xyz,1) , 'un',false ) )

% %visualise surfaces: endo, epi and myo lid
% cla
% patch('vertices',data(1).diastolic.endo.xyz,'faces',data(1).diastolic.endo.tri,'facecolor','none','EdgeColor','red')
% patch('vertices',data(1).diastolic.epi.xyz,'faces',data(1).diastolic.epi.tri,'facecolor','none','EdgeColor','blue')
% patch('vertices',myoB.xyz,'faces',myoB.tri,'facecolor','g')

%append epi surface, endo surface and myo lid
% Last_prevID  = size( EPI_ED.xyz , 1 );
% EPI_ED.xyz = [ EPI_ED.xyz ; B.xyz ];
% EPI_ED.tri = [ EPI_ED.tri ; B.tri + Last_prevID ];

dia_endo_prevID(i)  = size( transformed_data(1).diastolic.endo.xyz , 1 );
dia_epi_prevID(i) = size(myoB(i).dia.xyz,1);
myoB_full.xyz = [ transformed_data(1).diastolic.endo.xyz ; myoB.dia.xyz ; transformed_data(1).diastolic.epi.xyz ];
myoB_full.tri = [ transformed_data(1).diastolic.endo.tri; myoB.tri + dia_endo_prevID ; data(1).diastolic.epi.tri + dia_endo_prevID + dia_epi_prevID ];

%make sure that every triangle points outwards
myoB_full = FixNormals( myoB_full );

% visualise myocardium mesh
cla
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


end




