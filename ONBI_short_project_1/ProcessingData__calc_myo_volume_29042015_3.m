
 setenv('path',[getenv('path'),';','C:\Users\jesu2687\Documents\MATLAB\ErnestoCode\Tools\MESHES\vtk_libs']);
 addpath C:\Users\jesu2687\Documents\MATLAB\ErnestoCode\Tools\MESHES\
 addpath C:\Users\jesu2687\Documents\MATLAB\ErnestoCode\Tools\
 
cd C:\Users\jesu2687\Documents\MATLAB\ErnestoCode\

 % vtkCleanPolyData(EPI_ED) fix the possible replicated nodes and spurious
 % edges
for i = 1:400

    i
% transformed_data(i).diastolic.endo.xyz = reshape(transformed_data(i).diastolic.endo.xyz, size(data(i).diastolic.endo.xyz));
% transformed_data(i).diastolic.epi.xyz = reshape(transformed_data(i).diastolic.epi.xyz, size(data(i).diastolic.epi.xyz));
% transformed_data(i).systolic.endo.xyz = reshape(transformed_data(i).systolic.endo.xyz, size(data(i).systolic.endo.xyz));
% transformed_data(i).systolic.epi.xyz = reshape(transformed_data(i).systolic.epi.xyz, size(data(i).systolic.epi.xyz));

transformed_data(i).diastolic.endo.B = vtkFeatureEdges( vtkCleanPolyData(transformed_data(i).diastolic.endo) , 'BoundaryEdgesOn' , [] , 'FeatureEdgesOff' , [] );  %extracting the boundary
transformed_data(i).diastolic.endo.B.xyz = transformed_data(i).diastolic.endo.B.xyz( [2 1 3:end], : );  %fixing the connectivity.

transformed_data(i).diastolic.epi.B = vtkFeatureEdges( vtkCleanPolyData(transformed_data(i).diastolic.epi) , 'BoundaryEdgesOn' , [] , 'FeatureEdgesOff' , [] );  %extracting the boundary
transformed_data(i).diastolic.epi.B.xyz = transformed_data(i).diastolic.epi.B.xyz( [2 1 3:end], : );  %fixing the connectivity.

%load a manually produced myoB.tri
load('C:\Users\jesu2687\Documents\MATLAB\ONBI_short_project_1\myoB.mat')

%join the list of coordinates for endo and epi to be used to cover the myo
  myoB(i).dia.endo.xyz = transformed_data(i).diastolic.endo.B.xyz;
  myoB(i).dia.epi.xyz = transformed_data(i).diastolic.epi.B.xyz;
%   myoB(i).sys.endo.xyz = transformed_data(i).systolic.endo.B.xyz;
%   myoB(i).sys.epi.xyz = transformed_data(i).systolic.epi.B.xyz;
%   
  myoB(i).dia.xyz = [myoB(i).dia.endo.xyz ; myoB(i).dia.epi.xyz];
%   myoB(i).sys.xyz = [myoB(i).sys.endo.xyz ; myoB(i).sys.epi.xyz];

 
% %visualise surfaces: endo, epi and myo lid
% cla
% patch('vertices',data(1).diastolic.endo.xyz,'faces',data(1).diastolic.endo.tri,'facecolor','none','EdgeColor','red')
% patch('vertices',data(1).diastolic.epi.xyz,'faces',data(1).diastolic.epi.tri,'facecolor','none','EdgeColor','blue')
% patch('vertices',myoB.xyz,'faces',myoB.tri,'facecolor','g')

%append epi surface, endo surface and myo lid
% Last_prevID  = size( EPI_ED.xyz , 1 );
% EPI_ED.xyz = [ EPI_ED.xyz ; B.xyz ];
% EPI_ED.tri = [ EPI_ED.tri ; B.tri + Last_prevID ];

dia_endo_prevID  = 3267; %size( transformed_data(1).diastolic.endo.xyz , 1 );
dia_epi_prevID = 64;
total = dia_endo_prevID + dia_epi_prevID;

myoB_full(i).dia.xyz = [ transformed_data(i).diastolic.endo.xyz ; myoB(i).dia.xyz ; transformed_data(i).diastolic.epi.xyz ];
myoB_full(i).dia.tri = [ transformed_data(i).diastolic.endo.tri; myoB.tri + dia_endo_prevID ; transformed_data(i).diastolic.epi.tri + dia_endo_prevID + dia_epi_prevID ];

patch('vertices',myoB_full(i).dia.xyz,'faces',myoB_full(i).dia.tri,'facecolor','red')
%make sure that every triangle points outwards

% myoB_full(i).dia = FixNormals( myoB_full(i).dia );


% hold on
% plot3(myoB.xyz(:,1),myoB.xyz(:,2), myoB.xyz(:,3),'+')
% plot3(myoB_full.xyz(:,1),myoB_full.xyz(:,2), myoB_full.xyz(:,3),'o')

%calculate volume of myocardium
[myoVolumes(i).dia, myoCenterOfMasses(i).dia] = MeshVolume( myoB_full(i).dia );

myoDifference_volumes(i).dia = prod(diff( BBMesh( myoB_full(i).dia ) , 1  , 1 ) ) - myoVolumes(i).dia ;  %%it should be positive!!

end

%%
plot3(myoB(1).dia.xyz(:,1),myoB(1).dia.xyz(:,2), myoB(1).dia.xyz(:,3))
% visualise labelled points on the edge of the myocardium lid
plot3(myoB(1).dia.xyz(:,1) ,  myoB(1).dia.xyz(:,2) ,  myoB(1).dia.xyz(:,3) , '.r' )
text(myoB(1).dia.xyz(:,1) ,  myoB(1).dia.xyz(:,2) ,  myoB(1).dia.xyz(:,3) , arrayfun( @(id)sprintf('%d',id) , 1:size(myoB(1).dia.xyz,1) , 'un',false ) )


%% visualise myocardium mesh
cla
figure
patch('vertices',myoB_full(1).dia.xyz,'faces',myoB_full(1).dia.tri,'facecolor','red')

%% visualise center of mass within myocardium mesh
cla
patch('vertices',myoB_full(1).dia.xyz,'faces',myoB_full(1).dia.tri,'facecolor','g','facealpha',0.1);
hold on;
plot3( myoCenterOfMasses(1).dia(1) , myoCenterOfMasses(1).dia(2) , myoCenterOfMasses(1).dia(3) , '*r','markers',20 ); hold off




