%% Mesh utilities
% A mesh is a data structure to manage triangulations

%% Triangulation Data Structure
%
%  A triangulation is defined as a data structure of 2 matrices:
%  one to define the vertices coordinates:
%
%   *mesh.xyz*= [ x1  y1  z1  
%                 x2  y2  z2  
%                 x3  y3  z3  
%                 .
%                 .
%                 .
%                 xN  yN  zN  ]
% 
%  and the second define the triangles:
%
%   *mesh.tri* = [ 1 2 4 
%                  2 3 4
%                  .
%                  .
%                  vT1 vT2 vT3 ]  ( vertex 1 of triangle T )
% 
% 
% This means the triangle number 1 has vertices at points
%
%     (x1,y1,z1) ----- (x2,y2,z2)
%          |\              /| 
%           |\            /|
%            |\          /|
%             |\        /|
%              |\      /|
%               |\    /|
%              (x4,y4,z4)
% 
% 
%% An example
m.xyz=[ 0 0 0; 1 0 0; 2 1 0; 0 1 0; 1 1 0; 1 2 0];
m.tri=[ 1 2 4; 2 5 4; 2 3 5; 3 6 5];
plotMESH( m , 'numberpoints', 'numberelements')


%% Additional field
% It is very simple to add additional fields to the structure, for example fields which 
%   define scalar or vectorial data to the data. These fields could be defined for 
%   the vertices or for the triangles.
%   
% The fields which define properties over the nodes or vertices have to had the prefix
% 'xyz'. For example a field called 'xyzdata' will define values for each node.
% This have to had the data i the same order as the nodes.
%   
%   m.xyzdata= [ v1 ;  *the scalar value at the node 1*
%                v2 ;  *at node 2*
%                ...
%                vN ]  *at node N*
% 
%   
% Otherwise, the fields which names start with the prefix 'tri', define datas over the
% triangles. For example, 'trinormals':
%   
%   m.trinormals = [ nx1  ny1  nz1  ;  *define a vector for triangle 1*
%                    nx2  ny2  nz2  ;  *for triangle 2*
%                    ...
%                    nxT  nyT  nzT ]   *for triangle T*
%                  
                   
%% Landmarks 
% Landmarks are points defined for the meshes. The landmarks will be joined to the mesh
% for the transform operations.
%   
m.lmk  = [ 0.5 0.2 0 ; 1 1.5 0 ];
m.lmk2 = [ 2.0 1.5 0 ; 1.2 1.5 0 ];
plotMESH( m );
legend('show');


%% Transform Meshes
% Pose transformation operations
%   
%   sintax:
%     *trasformed_mesh = TransformMesh( original_mesh , pose_operations )*
%     *[trasformed_mesh, H] = TransformMesh( original_mesh , pose_operations )*
%     *[trasformed_mesh, H] = TransformMesh( ... , 'rotationUnits', ['Degrees' , 'Radians' ])*
%           the default is 'Degrees'
% 
%     *[trasformed_mesh, H] = TransformMesh( ... ,*
%                   *'Center', [ [xc yc zc] 'CenterMass' 'BoundingBox' 'Origin' ] )*
%           the default is 'Origin' ( (0,0,0) ). This is the point about
% the relative rotations ans scales will be applied
% 
%   
%   *pose_operations*: the concatenation of the follow basis operations
% 
%        *'TranslateX'*           , tx
%        *'TranslateY'*           , ty
%        *'TranslateZ'*           , tz
%        *'RotateX'*              , angle
%        *'RotateY'*              , angle
%        *'RotateZ'*              , angle
%        *'AxisAngle'*            , axis , angle       
%        *'RelativeRotateZ'*      , angle
%                 (rotate about Z axis about the defined center)
%        *'RelativeRotateX'*      , angle
%        *'RelativeRotateY'*      , angle
%        *'Scale'*                , factor
%        *'RelativeScale'*        , factor     
%                 (scale about the center of the defined center)
%        *'Pose3D'*               , [ tx,ty,tz,rry,rrx,rrz,s ]
%        *'CentertoCenterMass'*           
%                 (translate the centermass of the mesh to 0,0,0)
%        *'CentertoBoundingBox'*          
%                 (translate the center of the bounding box to 0,0,0)
%        *'TRANSFORM' ('H')*      , homogeneous_transform_matrix
% 
%
% An example:
%
plotMESH( m );
hold on
plotMESH( TransformMesh( m , 'tx',1,'ty',1,'rz',45,'s',2) , 'tc','r')
hold off

%% 
% The transformation are applied in the follow order, translate 1 on X , then , 
% translate 1 on Y, then rotate with respect to (0 0 0)
% and finally scale with respect to (0,0,0) by 2.
% 
% Another example:
%

plotMESH( m );
hold on
plotMESH( TransformMesh( m ,'center', 'bb','tx',1,'ty',1,'rrz',45,'s',2,'rz',45 ) , 'tc','r');
plotMESH( TransformMesh( m ,'center', 'bb','tx',1,'ty',1,'rrz',45,'rs',2,'rrz',45), 'tc','b');
hold off

%%
%  In this case some rotations and scale are done about the center of bounding box of the
%  original mesh.
% 
%  The landmarks are transformed together with the mesh.
% 

%% Plotting Meshes
%
%  plotmesh( mesh , 'Text'
%                   'NumberPoints'
%                   'NumberElements'
%                   'TriangleColor' , [r g b] or 'none' or LineSpec                 
%                   'EdgeColor'     , [r g b] or 'none' or LineSpec                 
%                   'PointData'     , data_values_on_points
%                   'TriangleData'  , data_values_on_triangles
% 
% The option 'Text' enforce the label of NumberPoints and NumberElements.
% 
% If PointData is a vectorial field of 3 component, then it will be draw as arrows. If it 
% is scalar data, then the surface will be painted. The same for the TriangleData.
% 
% Observation: the data over points has preference to be draw wrt the data over the triangles.
% 
% -
% 
% plMESH( mesh )
% 
% It is a interactive Viewer of VTK Polydata.
% 

%% Auxiliar Data Structures
% 
%   It is very straight forward to know, for each triangle, which points bellow to this.
%   This information are stored in the matrix tri, or connectivity matrix. The input 
%   for this matrix is the TRIANGLE_ID and the output is the POINTS_ID belonged to this.
%     
%   An usefull additional information for the triangulations, is to know, for each point,
%   which are the triangles who shared it. A function with the input, the POINT_ID and
%   an output with the TRIANGLE_ID is a inverse structure of the connectivity matrix.
%   
%   Unlike the number of nodes per element, which is a constant (3), the number of elements
%   surrounding a point can fluctuate from 1 in a corner to 8 or more in no smooth meshes.
% 
%   The most efficient way to store such varying data is using a linked list.
%   
%   The linked list is made up by two arrays, one store a list with the concatenation of
%   the TRIANGLE_IDs for each POINTS_ID and the second one is a list of pointers to the 
%   first one.
%   
%            pointer(p1)  pointer(p2) pointer(p3)        pointer(p4)
%            |            |           |                  |
%   list1: [ e1 e2 e3 e4  e2 e4 e5    e1 e3 e5 e6 e7 e8  .........]
%            \________ /  \______/    \_______________/  \____....
%                 |           |               |
%    POINTS_IDs   p1          p2              p3 ....
%    
%   These structure in the implementation is called ESUP (Elements SUrrounding Points )
% 
%   To create the structure, you have to use the follow command:
  
m= CreateESUP( m );
m.ESUP

%% 
plotMESH(m,'t')

%%
%   The pointer list is stored at m.ESUP.p and the element list at m.ESUP.el
%   
%   The element list for each point are:
%     m.ESUP.el( m.ESUP.p( POINT_ID )+1:m.ESUP.p( POINT_ID+1 ) )
%
% Examples:
elements_surrounding_point_2= m.ESUP.el( m.ESUP.p( 2 )+1:m.ESUP.p( 2+1 ) )

%%
%   This is implemented with the command ESUP
%   
%   sintax: elements_id= ESUP( points_id , mesh )
%   
%   If you specified more than one points_id, the command return the union of the ids
%   for each point_id.
%     
%   examples:
   
  ESUP( [1 3], m )
    
%%  
%   The same kind of structure is defined for the points. The command CreatePSUP, permit
%   to know the POINTS_IDs neighbours to one point. The command PSUP( points, mesh ) return
%   the POINT_IDs of the connected points to these.
%   
%   examples:
%   
  PSUP( 1 , m )
  
%%
  PSUP([1 2] , m )
  
%% Algorithms for Processing meshes
%   
%   *mesh_new= AppendMesh( mesh1, mesh2 )*
%   
%   Create in mesh_new a mesh with mesh1 and mesh2 appended. I think that the auxiliar 
%   fields and the attributes fields work well. If not, please tell me.
%   
% 
%   *convex_region_id= ConnectivityMesh( mesh )*
%   
%   Return a list of the convex_region_id for each triangle of the mesh.
%   
%   example:
m= CreateStructuredGrid(1:3,1:3);  
m= AppendMeshes( m , TransformMesh( m , 'ty', 2.5,'rrz', 45,'center','bb' ) );
m= AppendMeshes( m , TransformMesh( m , 'tx', 4.0,'rrz',-45,'center','bb' ) );
c= ConnectivityMesh( m );
plotMESH( m , 'td' , c ); colorbar

%%
%   *mesh_new= DeletePoints( points , mesh_original )*
%   
%   Remove the points *points* and the triangles surroungin these.
%   
%   example:
m= CreateStructuredGrid(1:3,1:3);  
plotMESH( m , 't');
hold on;
plotMESH( DeletePoints( [3 9], TransformMesh(m,'tx',2.5) ),'t');
hold off

%%
%   *mesh_new = MergePoints( mesh_original , threshold )*
%   
%   It merge the those points which are closer than *threshold*. Be carefull, you need
%   PLOC and SCAM to use it!!!

m= CreateStructuredGrid(1:3,1:3);
m= AppendMeshes( m , TransformMesh(m,'tx',2.1) );
plotMESH( m , 't'); axis off;
%%
plotMESH( MergePoints( CreatePLOC( CreateSCAM(m,1)),.15 ) , 't'); axis off;


%%
%   *mesh= MeshContour( X,Y,Z,C, value )*
%   
%   Return the isosurface mesh at the value *value* of the image *X,Y,Z,C*.
%   See isosurface in Matlab documentation.
%   
%   example:
  
[X,Y,Z]= meshgrid( -10:10 , -10:10 , -10:10 );
C= X*0;
C(:) = sqrt( X(:).^2 + abs( Y(:).^2.5 ) + 2*Z(:).^2 );
m= MeshContour( C,X,Y,Z, 15);
plotMESH( m , 'tc','b');
hold on;
plotMESH( MeshContour( C,X,Y,Z, 10) ,'tc','r' );
plotMESH( MeshContour( C,X,Y,Z, 20),'tc','g' );
hold off; hidden on
  
%%
%   *mesh_new= SubdivideMesh( mesh )*
%   
%   It divide each triangle into four. The new points will lay at the midle of
%   edges.
%   
%   example:
m= CreateStructuredGrid(1:2,1:2);
plotMESH( AppendMeshes( m , ...
            SubdivideMesh( TransformMesh( m ,'tx',1.5) ) ));

%%
%   *mesh_new= CleanMesh( mesh )*
%   
%   It merge the coincident points, remove the 0 area triangles and remove the unused points.
  
%%
%   *mesh_new= WarpMesh( mesh , origin_points, target_points, ...*
%                          *'Basis', 'R' or 'RLOGR'*  by default R
%                          *'Plot'*                   to visualize the warping
%   
%   It is an implementation of a TPS warp interpolation. The number of *origin_points*
%   has to be the same as the number of *target_points*.
%   For 3D data use the *R* basis which is the solution of the biharmonic eq. in 3D.
%   For 2D data, you could use *R*, but the solution of biharmonic is *RLOGR*.
%   
%   example:
m= CreateStructuredGrid(1:20,1:10);
m_new= WarpMesh( m , [ 0  0 0; 20 0 0; 0 10 0; 20 10 0;   5   5 0],...
                     [-1 -1 0; 20 1 0; 0 10 0; 24 14 0; 7.5 7.5 0],...
                     'plot','basis','r'); view(2)
%%
plotMESH( m , 'tc', 'none' );
hold on;
plotMESH( m_new , 'tc','r' );
hold off


%% 
% *mesh_new= SwapEdge( mesh , point1, point2 )*
% 
% Swap the edge which connect the points *point1* and *point2*.
% 
% example:
m     = CreateStructuredGrid(1:3,1:3);
plotMESH( m ,'t');
m_new = SwapEdge( m , 1,5 );
%%
plotMESH( m_new ,'t');

PSUP( 2 , m_new )
ESUP( 2 , m_new )   % the structures are automatically updated

%%
m_new = SwapEdge( m , 5,7 ); %no edge linking points 5 and 7... imposible to swap















