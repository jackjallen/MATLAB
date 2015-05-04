function DC = DiceCoefficient( A , B , varargin )



  volA = NaN;
  if iscell( A )
    volA = A{2};
    A    = A{1};
  end
  
  volB = NaN;
  if iscell( B )
    volB = B{2};
    B    = B{1};
  end
  
  [varargin,i,volA] = parseargs(varargin,'vola','$DEFS$',volA );
  [varargin,i,volB] = parseargs(varargin,'volb','$DEFS$',volB );
  
  [varargin,CLEAN] = parseargs(varargin,'noclean','$FORCE$',{false,true});
  

  if CLEAN
    A = struct( 'xyz' , A.xyz , 'tri' , A.tri );
    A = vtkCleanPolyData( A , 'SetAbsoluteTolerance',1e-10,'SetToleranceIsAbsolute',true );
    
    bounds = vtkFeatureEdges( A , 'BoundaryEdgesOn',[],'FeatureEdgesOff',[],'NonManifoldEdgesOff',[],'ManifoldEdgesOff',[]);
    if isfield( bounds , 'xyz' )
      warning('MESH A look open. try with vtkFillHolesFilter');
    end
  end
  
  if CLEAN
    B = struct( 'xyz' , B.xyz , 'tri' , B.tri );
    B = vtkCleanPolyData( B , 'SetAbsoluteTolerance',1e-10,'SetToleranceIsAbsolute',true );

    bounds = vtkFeatureEdges( B , 'BoundaryEdgesOn',[],'FeatureEdgesOff',[],'NonManifoldEdgesOff',[],'ManifoldEdgesOff',[]);
    if isfield( bounds , 'xyz' )
      warning('MESH B look open. try with vtkFillHolesFilter');
    end
  end
  

  I = MeshVolume( BooleanMeshes( A , 'intersection' , B , 'noclean' ) , 'noclean' );
  
  if isnan( volA ), volA = MeshVolume( A , 'noclean' ); end
  if isnan( volB ), volB = MeshVolume( B , 'noclean' ); end

  DC = 2*I/( volA + volB );

end
