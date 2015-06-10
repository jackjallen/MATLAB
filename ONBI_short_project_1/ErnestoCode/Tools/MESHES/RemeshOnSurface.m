function M = RemeshOnSurface( M , ISO , varargin )

  maxIT = 10;
  I     = [];
  elen  = 0;
  isovalue = 0;
  maxsmooth = Inf;
  

  
  [varargin,i,isovalue] = parseargs(varargin,'isoValue','$DEFS$', isovalue );
  [varargin,i,maxIT] = parseargs(varargin,'maxIT','$DEFS$', maxIT );
  [varargin,i,elen ] = parseargs(varargin,'EdgeLength','$DEFS$', elen );
  [varargin,i,I    ] = parseargs(varargin,'Image','$DEFS$', I );
  [varargin,i,resam] = parseargs(varargin,'Resample','$DEFS$', '-50e6' );  %-50e6
  [varargin,dibuja ] = parseargs(varargin,'plot','$FORCE$', 1 );
  [varargin,i,maxsmooth] = parseargs(varargin,'maxsmooth','$DEFS$', maxsmooth );
  %   dibuja = 1;


  if isa( ISO , 'I3D' )
    I = ISO;
    while 1
      try
        I.data = double( I.data );
        Vprintf('Extracting ISO ... ');
        if ischar( resam )
          Vprintf('( resample %s )', resam );
        else
          Vprintf('( resample %g )', resam );
        end
        ISO = isosurface( resample( I , resam ) ,  isovalue );
        Vprintf(' done\n');
        break;
      catch
        if ischar( resam )
          resam = num2str( str2double( resam ) / 10 );
        else
          resam = resam/10;
        end
      end
    end
  end

  
  M = FixMesh( M ); M = CleanMesh( M );
  ISO = FixMesh( ISO );
  vtkClosestElement( ISO );
  
  
  if dibuja
    figure;
    P = patch( 'vertices', M.xyz , 'faces' , M.tri , 'facecolor', [.6 .75 .75] , 'edgecolor', [ 0 0 0 ] );
    view(3);
    axis('equal');
    drawnow;
  end
  

  XYZp = M.xyz; XYZp(1) = XYZp(1) + eps( single( XYZp(1) ) );
  IT = 0;
  while IT < maxIT

    if isequal( XYZp , M.xyz )
      Vprintf('Converged\n');
      break;
    end
    IT = IT+1;
    Vprintf('IT: %4d of %d  ', IT , maxIT );
    
    
    XYZp = M.xyz;
    
    Vprintf('  (%d pts,%d tris)  ', size(M.xyz,1) , size(M.tri,1) );
    
    if elen > 0  % &&  ~rem( IT , 2 )
      M = CollapseSmallEdges( M , elen );
      Vprintf(' ->  (%d pts,%d tris)  ', size(M.xyz,1) , size(M.tri,1) );
    end

    if ~rem( IT , 5 )
      M = FixValences( M , 1000 );
    end

    XYZ = M.xyz;

    NORMALS = vtkComputeNormals( M , 'ComputeCellNormalsOff',[],'ConsistencyOn',[],'ComputePointNormalsOn',[],'AutoOrientNormalsOn',[]);
    
    %M = vtkWindowedSincPolyDataFilter(M,'SetNumberOfIterations',5,'SetFeatureAngle',180,'SetEdgeAngle',180,'FeatureEdgeSmoothingOn',[],'BoundarySmoothingOff',[],'GenerateErrorScalarsOff',[],'GenerateErrorVectorsOff',[],'NonManifoldSmoothingOn',[] );
    M = vtkSmoothPolyDataFilter(M,'SetNumberOfIterations',200,'SetFeatureAngle',180,'SetEdgeAngle',180,'FeatureEdgeSmoothingOn',[],'BoundarySmoothingOff',[],'GenerateErrorScalarsOff',[],'GenerateErrorVectorsOff',[]);
    
    
    D = M.xyz - XYZ;

    D = D - bsxfun(@times, sum( D.* NORMALS , 2 ) , NORMALS );
    
    DD = sqrt( sum( D.^2 , 2 ) );
    

    maxDD = max( DD(:) );
%     if maxDD > maxsmooth

      fprintf('\n max(DD) = %e  , disminuyendo smoothing  a %e\n', maxDD , maxsmooth );

      D( DD ~= 0 , : ) = bsxfun( @rdivide , D( DD ~= 0 , : ) , DD( DD ~= 0 ) );

      DD = min( DD , maxsmooth );
%       hist( DD , 100 ), drawnow

      M.xyz = XYZ + bsxfun( @times , D , DD );
%     end


    [kk,M.xyz] = vtkClosestElement( M.xyz );
    

    if ~isempty( I )
      Vprintf('  err: ( %g , %g ) ' , Mean_Var( I( M.xyz ) ) );
    end
    
    Vprintf('\n');


    if dibuja, set( P , 'vertices', M.xyz , 'faces' , M.tri ); pause(0.02); end
  
  end
  vtkClosestElement( 0 , 0 );
  
  function Vprintf( varargin )
    if maxIT > 1
      fprintf( varargin{:} );
    end
  end
  function mean_var = Mean_Var( x )
    mean_var = [ mean( x(:) ) ,  var( x(:) ) ];
  end

end

