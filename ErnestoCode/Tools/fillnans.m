function X = fillnans( X , Y )

  if ~isequal( size(X) , size(Y) )
    error('different sizes');
  end
  
  if ~strcmp( class(X) , class(Y) )
    error('different classes');
  end    

  w = ~isnan( X ) & ~isnan( Y );
  if ~isequal( X(w) , Y(w) )
    warning('warning!!! some non-NANS value are different!!. It will be avoided to overwrite the target!.');
  end
  
  
  w = isnan( X ) & ~isnan( Y );
  X(w) = Y(w);

end
