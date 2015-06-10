function M = BooleanMeshes( A , op , B , varargin )

  if nargin == 1
    M = [];
    conn = ConnectivityMesh( A );
    for c = unique( conn ).'
      MM = A;
      MM.tri( conn ~= c , : ) = [];
      
      if isempty( M )
        M = struct( 'xyz' , MM.xyz , 'tri' , MM.tri );
        M = vtkCleanPolyData( M , 'SetAbsoluteTolerance',1e-10,'SetToleranceIsAbsolute',true );

        M = checkFaces( M );
        
        bounds = vtkFeatureEdges( M , 'BoundaryEdgesOn',[],'FeatureEdgesOff',[],'NonManifoldEdgesOff',[],'ManifoldEdgesOff',[]);

        if isfield( bounds , 'xyz' )
          M = fillHoles( M );
        end
        
      else
        M = BooleanMeshes( M , 'union' , MM , 'clean','fill' );
      end
    end
    
    
    return;
  end



  [varargin,CLEAN] = parseargs(varargin,'noclean','$FORCE$',{false,true});
  [varargin,FILL ] = parseargs(varargin,'nofill' ,'$FORCE$',{false,true});
  if FILL, CLEAN = true; end

  switch lower(  op )
    case {'u','union','+'}
      op = 'UNION';
    case {'i','int','intersect','intersection'}
      op = 'INTERSECTION';
    case {'sd','symdiff','symmetricdifference'}
      op = 'SYMMETRIC_DIFFERENCE';
    case {'minus','-','a_minus_b'}
      op = 'A_MINUS_B';
    case {'b_minus_a'}
      op = 'B_MINUS_A';
    otherwise
      error('incorrect operation');
  end

  if CLEAN
    A = struct( 'xyz' , A.xyz , 'tri' , A.tri );
    A = vtkCleanPolyData( A , 'SetAbsoluteTolerance',1e-10,'SetToleranceIsAbsolute',true );

    B = struct( 'xyz' , B.xyz , 'tri' , B.tri );
    B = vtkCleanPolyData( B , 'SetAbsoluteTolerance',1e-10,'SetToleranceIsAbsolute',true );
  end

  
%   bounds = vtkFeatureEdges( A , 'BoundaryEdgesOn',[],'FeatureEdgesOff',[],'NonManifoldEdgesOff',[],'ManifoldEdgesOff',[]);
%   if isfield( bounds , 'xyz' ) && FILL
    [A,bounds] = fillHoles( A );
%   end
  if isfield( bounds , 'xyz' )
    warning('MESH A look open. try with vtkFillHolesFilter');
  end

  
%   bounds = vtkFeatureEdges( B , 'BoundaryEdgesOn',[],'FeatureEdgesOff',[],'NonManifoldEdgesOff',[],'ManifoldEdgesOff',[]);
%   if isfield( bounds , 'xyz' ) && FILL
    [B,bounds] = fillHoles( B );
%   end
  if isfield( bounds , 'xyz' )
    warning('MESH A look open. try with vtkFillHolesFilter');
  end

  
  if CLEAN
    A = checkFaces( A );
    B = checkFaces( B );
  end
  
  dirname = tmpname('BooleanMeshes_',8); 
  CLEANUP = onCleanup( @() rmdir(dirname,'s') );
  mkdir( dirname );

  A_fname = fullfile( dirname , 'A.obj' );
  B_fname = fullfile( dirname , 'B.obj' );
  M_fname = fullfile( dirname , 'M.vtk' );
  
  
  write_OBJ( A , A_fname );
  write_OBJ( B , B_fname );
  
  
  [p,f,e] = fileparts( mfilename('fullpath') );
  command = fullfile( p , 'carve' , 'intersect' );
  command = [ command , ' --vtk --triangulate --edge --rescale ' , '"' , strrep( A_fname , '\','\\') , '"' , '  ' , op , ' ' , '"' , strrep( B_fname , '\','\\') , '"' , ' > ' , '"' , M_fname , '"' ];
  
  [status,result] = system( command )
  if status
    error('no pudo!!!');
  else
    M = vtkPolyDataReader( M_fname );
  end



  function M = checkFaces( M )
    M.tri( any(M.tri == 0,2) | M.tri(:,1) == M.tri(:,2) | M.tri(:,1) == M.tri(:,3) | M.tri(:,2) == M.tri(:,3)  , : ) = [];
    
    
%     %flips = any( vtkPolyDataNormals(M,'SetComputePointNormals',false,'SetComputeCellNormals',true) ./ ComputeNormals( M ) < 0 , 2 );
%     flips = sum( vtkComputeNormals(M,'SetComputePointNormals',false,'SetComputeCellNormals',true) .* ComputeNormals( M ) , 2 ) < 0;
%     M.tri(flips,:) = M.tri(flips,[1 3 2]);
    M = vtkPolyDataNormals( M , 'ConsistencyOn',[],'AutoOrientNormalsOn',[],'ComputeCellNormalsOn',[],'ComputePointNormalsOff',[],'SplittingOff',[]);
  end
  
  function [M,bounds] = fillHoles( M )
%     Mp = M;
%     nM = vtkFillHolesFilter( M , 'SetHoleSize', 2 );
%     if isfield( nM , 'tri')
%       M = nM;
%     else
%       M = Mp;
%     end

    M = checkFaces( M );

    bounds = vtkFeatureEdges( M , 'BoundaryEdgesOn',[],'FeatureEdgesOff',[],'NonManifoldEdgesOff',[],'ManifoldEdgesOff',[]);

    if isfield( bounds , 'xyz' )
        %la malla tiene agujeros.
        error( 'error aca');
    end
    
  end    


end


