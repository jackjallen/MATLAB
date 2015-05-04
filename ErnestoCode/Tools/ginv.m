function [G,c] = ginv( X , c )

%   G = pinv(X);
%   return;

  if ndims(X)>2, error('a 2D matrix is expected'); end
  
  [n,k] = size(X);
  
  if n == k
    G = ( ( X \ eye(n) ) + ( eye(n) / X ) )/2;
    return;
  end
  
  TRANSPOSE = false;
  if n < k
    TRANSPOSE = true;
    X = X.';
  end
  
  if nargin > 1 &&  ischar(c)
    r = rank(X);
    [x,i] = sort(rand(1,n));
    [y,j] = sort(rand(1,k));

    C = X( i(1:r) , j(1:r) );

    while rank(C)<r
        [x,i] = sort(rand(1,n));
        [y,j] = sort(rand(1,k));
        C = X(i(1:r),j(1:r));
    end;

    C = C\eye(r); C = C';

    G = zeros(n,k);
    G(i(1:r),j(1:r)) = C;

    G = G.';
    
  else
    
    G = pinv( X );

    N = null(X.');
    if nargin < 2
      c = randn( size(N,2) , size(G,1) );
    end

    if ~isequal( c , 0 )
      G = G + ( N*c ).';
    end
      
  end
  
  if TRANSPOSE, G = G.'; end

end
