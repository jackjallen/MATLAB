function S = subspace_diff( P , R )


P = orth(P);
R = orth(R);

  N = null( P.'*R*R.' );
  S = orth( N*N.'*P );
              
end
