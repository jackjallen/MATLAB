function [ad,ADop] = adjoint( A )
% adjoint( A )*B(:) - vec( A*B - B*A )

  N = size(A,1);
  if ndims(A) ~= 2 || numel(A) ~= N*N, error('matrix must be square'); end
  
  ad = kron( eye(N)  , A ) - kron( A.' , eye(N) );
  
  if nargout > 1
    ADop = round( NumericalDiff( @(x)adjoint(x) , A ,'i' ) );
  end
end
