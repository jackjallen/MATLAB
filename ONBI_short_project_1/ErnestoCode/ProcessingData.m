
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




%%
 % vtkCleanPolyData(EPI_ED) fix the possible replicated nodes and spurious
 % edges
B = vtkFeatureEdges( vtkCleanPolyData(EPI_ED) , 'BoundaryEdgesOn' , [] , 'FeatureEdgesOff' , [] );  %extracting the boundary
B.xyz = B.xyz( [2 1 3:end], : );  %fixing the connectivity.

%center of the upper hole
C = mean( B.xyz , 1 );
B.xyz = [ B.xyz ; C ];

B.tri = [];
for t = 1:size(B.xyz,1)-2
  B.tri(end+1,:) = [ t , t+1 , size(B.xyz,1) ];
end
B.tri(end+1,:) = [ t+1 , 1 , size(B.xyz,1) ];

cla
patch('vertices',EPI_ED.xyz,'faces',EPI_ED.tri,'facecolor','r')
patch('vertices',B.xyz,'faces',B.tri,'facecolor','b')




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







