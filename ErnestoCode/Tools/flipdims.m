function x = flipdims( x , dims )

  if nargin < 2
    dims = find( size(x) > 1 );
  end

  for d = dims
    x = flipdim( x , d );
  end

end
