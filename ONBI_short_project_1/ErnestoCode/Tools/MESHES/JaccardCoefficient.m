function [ JC , dJC ] = JaccardCoefficient( A , B , varargin )

  
  if nargout <= 1

    %%mucho mas rapido porque no hay que calcular la union
    DC = DiceCoefficient( A , B , varargin{:} );
    JC = DC /( 2 - DC );
  
  else
    
    [ DC , dDC ] = DiceCoefficient( A , B , varargin{:} );
    JC = DC /( 2 - DC );
    dJC = 2 * dDC / ( 2 - DC )^2;
    
  end
  


%   [varargin,CLEAN] = parseargs(varargin,'noclean','$FORCE$',{false,true});
% 
% 
%   if CLEAN
%     A = struct( 'xyz' , A.xyz , 'tri' , A.tri );
%     A = vtkCleanPolyData( A , 'SetAbsoluteTolerance',1e-10,'SetToleranceIsAbsolute',true );
%     
%     bounds = vtkFeatureEdges( A , 'BoundaryEdgesOn',[],'FeatureEdgesOff',[],'NonManifoldEdgesOff',[],'ManifoldEdgesOff',[]);
%     if isfield( bounds , 'xyz' )
%       warning('MESH A look open. try with vtkFillHolesFilter');
%     end
%   end
%   
%   if CLEAN
%     B = struct( 'xyz' , B.xyz , 'tri' , B.tri );
%     B = vtkCleanPolyData( B , 'SetAbsoluteTolerance',1e-10,'SetToleranceIsAbsolute',true );
% 
%     bounds = vtkFeatureEdges( B , 'BoundaryEdgesOn',[],'FeatureEdgesOff',[],'NonManifoldEdgesOff',[],'ManifoldEdgesOff',[]);
%     if isfield( bounds , 'xyz' )
%       warning('MESH B look open. try with vtkFillHolesFilter');
%     end
%   end
%   
% 
%   I = MeshVolume( BooleanMeshes( A , 'intersection' , B , 'noclean' ) , 'noclean' );
%   U = MeshVolume( BooleanMeshes( A , 'union' ,        B , 'noclean' ) , 'noclean' );
% 
%   JC = I/U;


end
