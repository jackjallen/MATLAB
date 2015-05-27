
% loading a pair of meshes (for example from the subject id number 55 called SSM0119);
ENDO = load('Subject_55.mat','ENDO'); ENDO = ENDO.ENDO;
EPI  = load('Subject_55.mat','EPI' ); EPI  = EPI.EPI;

%by sure that .xyz and .tri fields are double class. It is the usual way
%that matlab manage the data. In case of using another datatype such as
%single or integer datatypes it is very probably to get a memory error and
%MATLAB crashes.
%In general, if you did not modify the datatype by hand in previous
%preprocess or in the reading steps, this is not necessary.
ENDO = struct( 'xyz' , double(ENDO.xyz) , 'tri' , double(ENDO.tri) );
EPI  = struct( 'xyz' , double(EPI.xyz)  , 'tri' , double(EPI.tri)  );


%compute the closest point from a 3D point X to the mesh ENDO
X = [ 0.1 , 0.2 , 0.3 ]   %be sure you haev 3 columns, each row specify a point coordinate
                          %again, X must be double, you can force it by:
                          % X = double( X )  although it is usually not
                          % necessary. If X is not double, unexpected
                          % result or MATLAB crashes will be obtained!
[ triID , xyz_closest_point , distance ] = vtkClosestElement( ENDO , X );
%"triID" is the triangle index of mesh ENDO where the closest point lies.
%"xyz_closest_point" are the coordinates of the point in the ENDO mesh (not
%necessarily a node) which is closest to X. If this point correspond to an
%edge or a node which are shared by several triangles, triID is one of them
%and I cannot ensure if it is the smaller or which of them.
%"distance" = sqrt( sum( (X - xyz_closest_point).^2 , 2 ) );

%X can contain several points and avoid a for loop on them
X = randn( 1000 , 3 )   %3 columns (x-, y- and z-coordinates), each row for each point
[ triIDs , xyz_closest_points , distances ] = vtkClosestElement( ENDO , X )


%Now, you can easily compute the distance from EPI 2 ENDO surfaces
[~,~,dEPI2ENDO] = vtkClosestElement( ENDO , EPI.xyz )

%and represent it as a scalar map on the EPI surface
figure; patch( 'vertices',EPI.xyz,'faces',EPI.tri,'facecolor','interp','cdata',dEPI2ENDO,'edgecolor',[1 1 1]*0.2)
axis equal;
view(3);
colormap jet
colorbar


%you can also subdivide the source points in order to get a finer
%resolution of the map. Each triangle will be divided into 4 by adding
%nodes at the midpoints of the edges
EPI_sub = SubdivideMesh( EPI , 2 ); %2 means that you are recursivelly subdividing the mesh 2 times.
figure; patch( 'vertices',EPI_sub.xyz,'faces',EPI_sub.tri,'facecolor','b','edgecolor',[1 1 1]*0.2);
axis equal;
view(3);

%or you can prefer to subdivide the mesh with the method documented in: http://www.vtk.org/doc/nightly/html/classvtkButterflySubdivisionFilter.html
EPI_sub = vtkButterflySubdivisionFilter( vtkButterflySubdivisionFilter( EPI ) );

figure; patch( 'vertices',EPI_sub.xyz,'faces',EPI_sub.tri,'facecolor','b','edgecolor',[1 1 1]*0.2);
axis equal;
view(3);

%and calculate the distance at a finer resolution
[~,~,dEPI2ENDO_sub] = vtkClosestElement( ENDO , EPI_sub.xyz );

%and represent it as a scalar map on the EPI surface
figure; patch( 'vertices',EPI_sub.xyz,'faces',EPI_sub.tri,'facecolor','interp','cdata',dEPI2ENDO_sub,'edgecolor','none')
axis equal;
view(3);
colormap jet
colorbar



%to do list!!
%unfold the EPI surface into a bull's eye plot (see for example:
%http://www.jcmr-online.com/content/pdf/1532-429X-10-S1-A229.pdf )
