function f = factor1(n)
% 
% factor1( 3999996938489837 )
% 

  if eps(n) > 1
    error('maximum precission reached');
  end

  f = 1;

  sn = ceil( sqrt(n) );
  
  P = [2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97];
  P( P > sn ) = [];

  ff = find( ~rem( n , P ) , 1 );
  if ~isempty(ff)
    
    f = P(ff);
    
  elseif sn > 100
    
    P = primes( sn );
    P = P( P > 100 );

    ff = find( ~rem( n , P ) , 1 );
    if ~isempty(ff)
      f = P(ff);
    end
    
  end
  
end
