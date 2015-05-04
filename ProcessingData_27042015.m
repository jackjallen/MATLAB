 setenv('path',[getenv('path'),';','C:\Users\jesu2687\Documents\MATLAB\ErnestoCode\Tools\MESHES\vtk_libs']);
 addpath C:\Users\jesu2687\Documents\MATLAB\ErnestoCode\Tools\MESHES\
 addpath C:\Users\jesu2687\Documents\MATLAB\ErnestoCode\Tools\
 
cd C:\Users\jesu2687\Documents\MATLAB\ErnestoCode\

%reading the vertices coordinates
fid = fopen( 'Subject1\SSM0001.ED.endo.vertices.csv','r');
EPI_ED.xyz = textscan( fid , '%f,%f,%f' );
EPI_ED.xyz = cell2mat( EPI_ED.xyz );
fclose( fid );

%reading triangles
fid = fopen( 'Subject1\TriangleFaces.csv ','r');
EPI_ED.tri = textscan( fid , '%f,%f,%f' );
EPI_ED.tri = cell2mat( EPI_ED.tri );
fclose( fid );
%%


load('C:\Users\jesu2687\Documents\MATLAB\ONBI_short_project_1\project1_data.mat');
%%
 % vtkCleanPolyData(EPI_ED) (JA:function provided by Ernesto)fix the possible replicated nodes and spurious
 % edges
for i = 1:400 %calculate the hole edge line coordinates xyz B for all patients
%diastolic, endo
data(i).diastolic.endo.B = vtkFeatureEdges( vtkCleanPolyData(data(i).diastolic.endo) , 'BoundaryEdgesOn' , [] , 'FeatureEdgesOff' , [] );  %extracting the boundary
data(i).diastolic.endo.B.xyz = data(i).diastolic.endo.B.xyz( [2 1 3:end], : );  %fixing the connectivity.
%diastolic, epi
data(i).diastolic.epi.B = vtkFeatureEdges( vtkCleanPolyData(data(i).diastolic.epi) , 'BoundaryEdgesOn' , [] , 'FeatureEdgesOff' , [] );  %extracting the boundary
data(i).diastolic.epi.B.xyz = data(i).diastolic.epi.B.xyz( [2 1 3:end], : );  %fixing the connectivity.
%systolic, endo
data(i).systolic.endo.B = vtkFeatureEdges( vtkCleanPolyData(data(i).systolic.endo) , 'BoundaryEdgesOn' , [] , 'FeatureEdgesOff' , [] );  %extracting the boundary
data(i).systolic.endo.B.xyz = data(i).systolic.endo.B.xyz( [2 1 3:end], : );  %fixing the connectivity.
%systolic, epi
data(i).systolic.epi.B = vtkFeatureEdges( vtkCleanPolyData(data(i).systolic.epi) , 'BoundaryEdgesOn' , [] , 'FeatureEdgesOff' , [] );  %extracting the boundary
data(i).systolic.epi.B.xyz = data(i).systolic.epi.B.xyz( [2 1 3:end], : );  %fixing the connectivity.

%find the center of the upper holes (e.g diastolic_endo_C)

% diastolic, endo
diastolic_endo_C = mean( data(i).diastolic.endo.B.xyz , 1 );
data(i).diastolic.endo.B.xyz = [ data(i).diastolic.endo.B.xyz ; diastolic_endo_C ];
% diastolic, epi
diastolic_epi_C = mean( data(i).diastolic.epi.B.xyz , 1 );
data(i).diastolic.epi.B.xyz = [ data(i).diastolic.epi.B.xyz ; diastolic_epi_C ];
% systolic, endo
systolic_endo_C = mean( data(i).systolic.endo.B.xyz , 1 );
data(i).systolic.endo.B.xyz = [ data(i).systolic.endo.B.xyz ; systolic_endo_C ];
% systolic, epi
systolic_epi_C = mean( data(i).systolic.epi.B.xyz , 1 );
data(i).systolic.epi.B.xyz = [ data(i).systolic.epi.B.xyz ; systolic_epi_C ];



    data(i).diastolic.endo.B.tri = [];
    data(i).diastolic.epi.B.tri = [];
    data(i).systolic.endo.B.tri = [];
    data(i).systolic.epi.B.tri = [];
    
    %diastolic, endo
    for t = 1:size(data(i).diastolic.endo.B.xyz,1)-2
        data(i).diastolic.endo.B.tri(end+1,:) = [ t , t+1 , size(data(i).diastolic.endo.B.xyz,1) ];
    end
   data(i).diastolic.endo.B.tri(end+1,:) = [ t+1 , 1 , size(data(i).diastolic.endo.B.xyz,1) ];
    
    %diastolic, epi
    for t = 1:size(data(i).diastolic.epi.B.xyz,1)-2
       data(i).diastolic.epi.B.tri(end+1,:) = [ t , t+1 , size(data(i).diastolic.epi.B.xyz,1) ];
    end
    data(i).diastolic.epi.B.tri(end+1,:) = [ t+1 , 1 , size(data(i).diastolic.epi.B.xyz,1) ];
    
     %systolic, endo
    for t = 1:size(data(i).systolic.endo.B.xyz,1)-2
       data(i).systolic.endo.B.tri(end+1,:) = [ t , t+1 , size(data(i).systolic.endo.B.xyz,1) ];
    end
    data(i).systolic.endo.B.tri(end+1,:) = [ t+1 , 1 , size(data(i).systolic.endo.B.xyz,1) ];
    
     %systolic, epi
    for t = 1:size(data(i).systolic.epi.B.xyz,1)-2
       data(i).systolic.epi.B.tri(end+1,:) = [ t , t+1 , size(data(i).systolic.epi.B.xyz,1) ];
    end
    data(i).systolic.epi.B.tri(end+1,:) = [ t+1 , 1 , size(data(i).systolic.epi.B.xyz,1) ];

%Plot the lid and the body of the LV
cla
patch('vertices',data(1).systolic.epi.xyz,'faces',data(1).systolic.epi.tri,'facecolor','r')
patch('vertices',data(1).systolic.epi.B.xyz,'faces',data(1).systolic.epi.B.tri,'facecolor','b')

%append both meshes
sys_epi_Last_prevID  = size( data(i).systolic.epi.xyz , 1 );
sys_endo_Last_prevID  = size( data(i).systolic.endo.xyz , 1 );
dia_epi_Last_prevID  = size( data(i).diastolic.epi.xyz , 1 );
dia_endo_Last_prevID  = size( data(i).diastolic.endo.xyz , 1 );

data(i).systolic.epi.xyz = [ data(i).systolic.epi.xyz ; data(i).systolic.epi.B.xyz ];
data(i).systolic.endo.xyz = [ data(i).systolic.endo.xyz ; data(i).systolic.endo.B.xyz ];
data(i).diastolic.epi.xyz = [ data(i).diastolic.epi.xyz ; data(i).diastolic.epi.B.xyz ];
data(i).diastolic.endo.xyz = [ data(i).diastolic.endo.xyz ; data(i).diastolic.endo.B.xyz ];

data(i).systolic.epi.tri = [ data(i).systolic.epi.tri ; data(i).systolic.epi.B.tri + Last_prevID ];
data(i).systolic.endo.tri = [ data(i).systolic.endo.tri ; data(i).systolic.endo.B.tri + Last_prevID ];
data(i).diastolic.epi.tri = [ data(i).diastolic.epi.tri ; data(i).diastolic.epi.B.tri + Last_prevID ];
data(i).diastolic.endo.tri = [ data(i).diastolic.endo.tri ; data(i).diastolic.endo.B.tri + Last_prevID ];

%make sure that every triangle points outwards
data(i).systolic.epi = FixNormals( data(i).systolic.epi );
data(i).systolic.endo = FixNormals( data(i).systolic.endo );
data(i).diastolic.epi = FixNormals( data(i).diastolic.epi );
data(i).diastolic.endo = FixNormals( data(i).diastolic.endo );

%patch('vertices',EPI_ED.xyz,'faces',EPI_ED.tri,'facecolor','r')

%calculate LV volume and center of mass (CoM)
[data(i).systolic.epi.volume, data(i).systolic.epi.centerofmass] = MeshVolume(data(i).systolic.epi);
[data(i).systolic.endo.volume, data(i).systolic.endo.centerofmass] = MeshVolume(data(i).systolic.endo);
[data(i).diastolic.epi.volume, data(i).diastolic.epi.centerofmass] = MeshVolume(data(i).diastolic.epi);
[data(i).diastolic.endo.volume, data(i).diastolic.endo.centerofmass] = MeshVolume(data(i).diastolic.endo);

data(i).systolic.epi.differencevolume = prod(diff( BBMesh( EPI_ED ) , 1  , 1 ) ) - Volume;  %%it shoud be positive!!
data(i).systolic.endo.differencevolume = prod(diff( BBMesh( EPI_ED ) , 1  , 1 ) ) - Volume;  %%it shoud be positive!!
data(i).diastolic.epi.differencevolume = prod(diff( BBMesh( EPI_ED ) , 1  , 1 ) ) - Volume; %%it shoud be positive!!
data(i).diastolic.endo.differencevolume = prod(diff( BBMesh( EPI_ED ) , 1  , 1 ) ) - Volume; %%it shoud be positive!!


%cla
%patch('vertices',EPI_ED.xyz,'faces',EPI_ED.tri,'facecolor','b','facealpha',0.2); hold on; plot3( CenterOfMass(1) , CenterOfMass(2) , CenterOfMass(3) , '*r' ); hold off

%%

%plot3( B.xyz(:,1) ,  B.xyz(:,2) ,  B.xyz(:,3) , '.r' )
%text(  B.xyz(:,1) ,  B.xyz(:,2) ,  B.xyz(:,3) , arrayfun( @(id)sprintf('%d',id) , 1:size(B.xyz,1) , 'un',false ) )

%% JA: plotting the volumees...

   sys_epi_volumes(i,1) = data(i).systolic.epi.volume
   sys_endo_volumes(i,1) = data(i).systolic.endo.volume
   dia_epi_volumes(i,1) = data(i).diastolic.epi.volume
   dia_endo_volumes(i,1) = data(i).diastolic.endo.volume

end
   
hold on
plot(sys_epi_volumes)
plot(sys_endo_volumes)
plot(dia_epi_volumes)
plot(dia_endo_volumes)
legend('sys epi volumes', 'sys endo volumes','dia epi volumes','dia endo volumes')





