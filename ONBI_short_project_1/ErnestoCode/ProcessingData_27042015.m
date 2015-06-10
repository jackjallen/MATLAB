
setenv('path',[getenv('path'),';','F:\ErnestoCode\Tools\MESHES\vtk_libs']);
addpath F:\ErnestoCode\Tools\MESHES\
addpath F:\ErnestoCode\Tools\

cd F:\ErnestoCode\

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
clear all
close all
clc



%%
 % vtkCleanPolyData(EPI_ED) (JA:function provided by Ernesto)fix the possible replicated nodes and spurious
 % edges
i = 1 %for i = 1:400 %for all patients
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

%end

%center of the upper holes
for i = 1:400 % all patients
% diastolic, endo
diastolic_endo_C = mean( data(i).diastolic.endo.B.xyz , 1 );
data(i).diastolic.endo.B.xyz = [ data(i).diastolic.endo.B.xyz ; diastolic_endo_C ];
% diastolic, epi
diastolic_epi_C = mean( data(i).diastolic.epi.B.xyz , 1 );
data(i).diastolic.epi.B.xyz = [ data(i).diastolic.epi.B.xyz ; diastolic_epi_C ];
% systolic, endo
systolic_endo_C = mean( data(i).diastolic.epi.B.xyz , 1 );
data(i).diastolic.epi.B.xyz = [ data(i).diastolic.epi.B.xyz ; diastolic_epi_C ];
% systolic, epi
systolic_epi_C = mean( data(i).diastolic.epi.B.xyz , 1 );
data(i).diastolic.epi.B.xyz = [ data(i).diastolic.epi.B.xyz ; diastolic_epi_C ];
end


for i = 1:400
    
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

end


cla
patch('vertices',data(1).systolic.epi.xyz,'faces',data(1).systolic.epi.tri,'facecolor','r')
patch('vertices',data(1).systolic.epi.B.xyz,'faces',data(1).systolic.epi.B.tri,'facecolor','b')




%append both meshes
Last_prevID  = size( EPI_ED.xyz , 1 );
EPI_ED.xyz = [ EPI_ED.xyz ; B.xyz ];
EPI_ED.tri = [ EPI_ED.tri ; B.tri + Last_prevID ];

%make sure that every triangle points outwards
EPI_ED = FixNormals( EPI_ED );

patch('vertices',EPI_ED.xyz,'faces',EPI_ED.tri,'facecolor','r')


[Volume,CenterOfMass] = MeshVolume( EPI_ED )

difference_volume = prod(diff( BBMesh( EPI_ED ) , 1  , 1 ) ) - Volume   %%it shoud be positive!!

cla
patch('vertices',EPI_ED.xyz,'faces',EPI_ED.tri,'facecolor','b','facealpha',0.2); hold on; plot3( CenterOfMass(1) , CenterOfMass(2) , CenterOfMass(3) , '*r' ); hold off

%%

plot3( B.xyz(:,1) ,  B.xyz(:,2) ,  B.xyz(:,3) , '.r' )
text(  B.xyz(:,1) ,  B.xyz(:,2) ,  B.xyz(:,3) , arrayfun( @(id)sprintf('%d',id) , 1:size(B.xyz,1) , 'un',false ) )







