function I = Mesh2I3D( M , I , varargin )

  [varargin,DIST] = parseargs( varargin , 'dist' , '$FORCE$' , {true,false} );

  
  M = FixMesh( M );

  
  if isa( I , 'I3D' )  &&  ~DIST
    I.data = false( size(I,1:3) );
    
    M.xyz = transform( M.xyz , maketransform( I.SpatialTransform ,'inv') );

    sz = [  numel( find( I.X > min( M.xyz(:,1) ) , 1 ,'first' ):find( I.X < max( M.xyz(:,1) ) , 1 ,'last' ) )  ...
            numel( find( I.Y > min( M.xyz(:,2) ) , 1 ,'first' ):find( I.Y < max( M.xyz(:,2) ) , 1 ,'last' ) )  ...
            numel( find( I.Z > min( M.xyz(:,3) ) , 1 ,'first' ):find( I.Z < max( M.xyz(:,3) ) , 1 ,'last' ) )  ];
    [sz,id] = min( sz );
    switch id
      case 1
        [Y,Z] = ndgrid( I.Y , I.Z );
        for i = find( I.X > min( M.xyz(:,1) ) , 1 ,'first' ):find( I.X < max( M.xyz(:,1) ) , 1 ,'last' )
          [slice,conn] = SliceMesh( M , [ I.X(i) 0 0  ; 1 0 0 ] , true );
          I.data(i,:,:) = reshape( inpoly( [ Y(:).' ; Z(:).' ] , [ slice(:,2).' ; slice(:,3).' ] , conn.' ), size(I,[2 3]) );
        end
        
      case 2
        [X,Z] = ndgrid( I.X , I.Z );
        for j = find( I.Y > min( M.xyz(:,2) ) , 1 ,'first' ):find( I.Y < max( M.xyz(:,2) ) , 1 ,'last' )
          [slice,conn] = SliceMesh( M , [ 0 I.Y(j) 0  ; 0 1 0 ] , true );
          I.data(:,j,:) = reshape( inpoly( [ X(:).' ; Z(:).' ] , [ slice(:,1).' ; slice(:,3).' ] , conn.' ), size(I,[1 3]) );
        end
        
      case 3
        XY = ndmat( I.X , I.Y ).';
        for k = find( I.Z > min( M.xyz(:,3) ) , 1 ,'first' ):find( I.Z < max( M.xyz(:,3) ) , 1 ,'last' )
          [slice,conn] = SliceMesh( M , [ 0 0 I.Z(k) ; 0 0 1 ] , true );
          if isempty(slice), continue; end
          I.data(:,:,k) = reshape( inpoly( XY , [ slice(:,1).' ; slice(:,2).' ] , conn.' ), size(I,[1 2]) );
        end
    end

  elseif isa( I , 'I3D' )  &&  DIST
    
    [kk,kk,I.data(:)] = vtkClosestElement( M , I.XYZ );
    
    interior = mesh2I3D( M , I );
    I.data( interior.data ) = - I.data( interior.data );
    
  elseif isscalar( I ) || ( isnumeric( I ) && numel( I ) == 3 )
    
    if numel( I ) == 1, I = [ I I I]; end
    
    bb = [ min( M.xyz , [] , 1 ) ; max( M.xyz , [] , 1 ) ];
    bb = [ bb(1,:) - eps( bb(1,:) )*100 ; bb(2,:) + eps( bb(2,:) )*100 ];
    
    I = I3D( 'X' , linspace( bb(1,1) , bb(2,1) , I(1) )  ,...
             'Y' , linspace( bb(1,2) , bb(2,2) , I(2) )  ,...
             'Z' , linspace( bb(1,3) , bb(2,3) , I(3) )  );
    I = padding( I , [1 1 1] , 'value' , 0 );
    
    if DIST
      I = mesh2I3D( M , I , 'dist' );
    else
      I = mesh2I3D( M , I );
    end
    
  end


end
