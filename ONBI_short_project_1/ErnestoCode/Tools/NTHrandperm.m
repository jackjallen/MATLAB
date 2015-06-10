function p = NTHrandperm( nth , n , fast )

  if nargin < 3
    fast = false;
  elseif strcmpi( fast  , 'fast' )
    fast = true;
  elseif strcmpi( fast  , 'slow' )
    fast = false;
  else
    error('third argument must be ''fast'' or ''slow''');
  end
    
    

  if nth == 1
    p = 1:n;
    return;
  end

  if fast
    
    rand( 'twister' , nth*(n+1) );
    [ignore,p] = sort( rand(1,n) );
    
  else
  
    s_old = [];
    try
      s_old = randseed;
    end

    rand( 'twister' , nth*(n+1) );
    [ignore,p] = sort( rand(1,n) );

    if ~isempty( s_old )
      rand( s_old{:} );
    end
    
  end
  
end
