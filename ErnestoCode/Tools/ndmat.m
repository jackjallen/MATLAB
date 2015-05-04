function coords = ndmat( varargin )
%
% create a n-dimensional grid
%
% gridd( x0:x1 , y0:y1 , z0:z1 , w0:w1 , ..., ['Sort'] )
%
% gridd( [dimx dimy ... dimz] )
%   goes from 1:dimx , 1:dimy ....
%
%

  [varargin,ordenar]= parseargs( varargin,'Sort', '$FORCE$',1 );
  [varargin,nocat  ]= parseargs( varargin,'nocat',    '$FORCE$',1 );
  [varargin,ascells]= parseargs( varargin,'asCELL','asCELLS',    '$FORCE$',1 );
  [varargin,i,comp ]= parseargs( varargin,'Comp' );

  if iscell( varargin{1} )
    varargin = cellfun( @(x) 1:x , varargin{1} , 'UniformOutput',false );
  end
  ndims= numel( varargin );


  if ordenar
    varargin = cellfun( @(x) sort( x(:) ) , varargin , 'UniformOutput',false );
  else
    varargin = cellfun( @(x)       x(:)   , varargin , 'UniformOutput',false );
  end
  
  if ndims == 1
    coords = varargin;
  else
    [ coords{1:ndims} ]= ndgrid( varargin{:} );
  end


  if ~isempty(comp)
    coords = coords{comp};
    return;
  end


  if nocat && ~ascells
    coords = cat( ndims+1 , coords{:} );
  elseif nocat && ascells
    
  elseif ~nocat && ~ascells
    coords = cellfun( @(x) x(:) , coords ,'UniformOutput',false );
    coords = cat( 2 , coords{:} );
  elseif ~nocat && ascells
    coords = cellfun( @(x) x(:) , coords ,'UniformOutput',false );
    coords = cat( 2 , coords{:} );
    coords = mat2cell( coords , ones( size(coords,1) , 1) , size(coords,2) );
  end

  
  
end
