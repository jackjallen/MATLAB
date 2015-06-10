function x = shuffle( x , fast )
  
  if nargout < 2
    fast = false;
  end

  n = numel(x);
  if n == 1,
    n = prod(size(x));
  end
  
  if fast

    x(:) = x( randperm( n ) );
    
  else
  
    s_old = randseed;

    x(:) = x( randperm( n ) );

    rand( s_old{:} );
  end

end
