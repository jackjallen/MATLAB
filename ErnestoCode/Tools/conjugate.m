function Q = conjugate( Q , C , useInverse )

  if nargin < 3 , useInverse = false; end
  if ~isequal( size( Q ) , size( C ) ), error('size must be equal'); end

  if ~useInverse
    
    Q = ( C * Q ) / C;
    
  else
    
    Q = ( C \ Q ) * C;
    
  end

end
