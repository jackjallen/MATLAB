function h = estimate_pdf( x , N )
%{

N = 500;                        %Number of samples
x = [0.5-.5*log(rand(1,N))...   %Complex
    0.4 + 0.03*randn(1,N/5) ...
    1.2 + 0.03*randn(1,N/10)];

hh = estimate_pdf( x , 1000 );
plot( hh(:,1) , hh(:,2) ,'r' );

%}


  if ~isnumeric(x), error('x must be a numerical array'); end
  if  numel(x) < 3, error('x have at least 3 elements'); end
  if  ~isvector(x), error('x must be a vector'); end

  if nargin < 2 || isempty(N)
    N=1000;
  end
  
  if isscalar(N)
    m = min(x(:));
    M = max(x(:));
    r = (M+m)/2 + 1.5*[-1 1]*(M-m)/2;
    
    N = linspace( r(1) , r(2) , N );
  end
  if ~isvector(N) , error('N must be an increasing vector or a scalar'); end
  if any(diff(N)<0) , error('N must be an increasing vector or a scalar'); end
  

  [h(:,2),h(:,1)] = ssvkernel(x(:),N(:).');

end

