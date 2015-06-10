function result = is3dtform( M , type )

  if nargin == 0, error('at least 1 argument is required'); end

  tol = 0;

  result = false;
  if ~isnumeric( M )
    return;
  end
  if ndims( M ) ~= 2
    return; 
  end
  if size( M , 1 ) ~= 4
    return; 
  end
  if size( M , 2 ) ~= 4
    return; 
  end
  if abs( M(4,4) - 1 ) > tol
    return; 
  end
  if any( abs( M(4,1:3) ) > tol )
    return; 
  end

  if nargin < 2
    type = 'affine';
  end


  tol = tol + 10*eps(1);
  switch lower( type )
    case {'affine'}
      
    case {'rigid'}
      M = logm( M(1:3,1:3) );
      if any( abs( M([1 5 9]) ) > tol )
        return;
      end
      if any( abs( M([2 3 6]) + M([4 7 8]) ) > tol )
        return; 
      end
      
    case {'sim','similarity'}
      M = logm( M(1:3,1:3) );
      if any( abs( M([2 3 6]) + M([4 7 8]) ) > tol )
        return; 
      end
      
    case {'traslation','translation','t'}
      if any( abs( M(1:3,1:3) - eye(3) ) > tol )
        return;
      end
      
    otherwise
      error('unrecognizable type');
  end
  
  
  result = true;
end
