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

  [varargin,ordenar]= parseargs( varargin,'sort','s', '$FORCE$',{1,0} );
  [varargin,nocat  ]= parseargs( varargin,'nocat',    '$FORCE$',{1,0} );

  if iscell( varargin{1} )
    varargin= varargin{:};
  end
  ndims= numel( varargin );

  if ndims == 1
    dims = varargin{1};
    ndims = numel( dims );
    for d=1:ndims
      varargin{d}= 1:dims(d);
    end
  end

  if ordenar
    for d=1:ndims
      varargin{d}= sort( varargin{d} );
    end
  end

  [ coords{1:ndims} ]= ndgrid( varargin{:} );

  if nocat
    coords = cat( ndims+1 , coords{:} );
  else
    coords = cellfun( @(x) x(:) , coords ,'UniformOutput',false );
    coords = cat( 2 , coords{:} );
  end

end
