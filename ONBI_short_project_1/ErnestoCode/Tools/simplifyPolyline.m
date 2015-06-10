function [xy,yy] = simplifyPolyline( xy , tol )

  if size( xy , 1 ) < 3, return; end
  if nargin < 2, tol = 1e-10; end

  
  m = diff( xy , 1 , 1 );
  m = noinfs( m , NaN );
  m = atan2( m(:,1) , m(:,2) );
  
  m = ~~[ 1 ; ~( abs( diff( m , 1 , 1 ) ) <= tol ) ; 1 ];
  
  
  xy = xy( m , : );
  
  if nargout > 1
    yy = xy(:,2);
    xy = xy(:,1);
  end

end
