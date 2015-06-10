function x = crandn( n , m , s )
% 
%   x = crandn( n , MU , SIGMA )
%   
%   return N  row-vectors with mean MU and covariance SIGMA
% 


  if ~isscalar(n) || n < 0, error('N has to be a positive scalar'); end
  if nargin < 2, D = 1; m = 0; s = 1; end
  if isempty(m)
    if nargin < 3
      D = 1;
    else
      D = size(s,1);
    end
    m = zeros([1,D]);
  else
    D = numel(m);
  end
  if nargin < 3
    s = eye(D);
  end
  
  
  if size(m,1) ~= 1, error('M is expected as a row'); end
  if ndims(s) ~= 2 || size(s,1) ~= size(s,2) || maxnorm( s - s.' ) > 1e-8
    error('S is expected as a square symmetric matrix');
  end
  if numel(m) ~= size(s,1)
    error('M and S have to have the same dimension!!');
  end

  if n == 0, x = zeros([0,D]); return; end

  x = randn(n,D);

  C = chol( s );
  x = bsxfun( @plus , x*C , m );

end