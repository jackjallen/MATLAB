function [coeffs, res_var] = genGLM(y, M , G)
% [coeffs, res_var] = genGLM(y, M , G)

sz_y = size(y);

G0 = null(G);
[Q,R] = qr( M*G0 , 0 ); 
coeffs = G0 * linsolve( R , Q.'*y(:,:) , struct('UT',true) ); 


if nargout>1
  if isa(y, 'double')
    res_var = reshape(residualVariance(y(:,:), M, coeffs), [sz_y(2:end) 1]); 
  else
    res_var = reshape(residualVariance_single(y(:,:), single(M), single(coeffs)), [sz_y(2:end) 1]); 
  end
end

coeffs = reshape(coeffs, [size(G,2) sz_y(2:end)]);