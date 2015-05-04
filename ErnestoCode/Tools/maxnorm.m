function a_ = maxnorm( a , b )

  if nargin > 1, a = a - b; end

  a= max( abs( a(:) ) );
  
  if nargout > 0
    a_ = a;
  else
    fprintf('%.15g\n',a);
  end

end
