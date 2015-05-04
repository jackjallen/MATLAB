function mpiX = mpinv( X )
% MPINV Moore-Penrose inverse matrix (Full matrices)
%
% Reference:
%   Katsikis, Pappas, Petralias, 
%   An improved method for the computation of the Moore Penrose inverse matrix, 
%   Applied Mathematics and Computation, 217, (2011) pp. 9828-9834.

  m = size(X,1);
  [Q,R,P] = qr(X);
  r = sum( any( abs(R) > 1e-5 , 2 ) );
  
  
  if m > r,

    R = R(1:r,:).';
    
    C = R * R.';
    R = C \ R;

    Q = Q(:,1:r);
    
  else
    
    C = R * R.';
    R = C \ R;
    R = R.';
    
  end
  

  mpiX = ( P *  R ) * Q.';
  
end
