function x = getv( x , varargin )

  if numel( varargin ) == 1 && iscell( varargin{1} )
    
    x = x{ varargin{1}{1} };

  elseif isstruct( x ) && strcmp( varargin{1} , '.' )
    
    x = x.(varargin{2});
    
  elseif isstruct( x ) && ~strcmp( varargin{1} , '.' ) && strncmp( varargin{1} , '.' , 1 )

    x = x.(varargin{1}(2:end));
    
  else

    x = x( varargin{:} );

  end

end
