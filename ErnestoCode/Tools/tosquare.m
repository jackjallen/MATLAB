function x = tosquare( x )

  P = numel(x);
  N = sqrt(P);
  
  if mod(N,1), error('not squared!!'); end
  
  x = reshape( x , [ N , N ] );

end
