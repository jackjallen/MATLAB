function L = logmsym( M )

  if isa(M,'sym')
    [v,d] = eig(M);
    L = v * diag( log( diag(d) ) ) / v;
  else
    L = logm(M);
  end
    
end
