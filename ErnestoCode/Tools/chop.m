function x = chop( x , thres )

  if nargin < 2
    thres = 10*eps(1);
  end
  
  if thres > 0

    x( abs( x ) <= thres ) = 0;
    
  elseif thres < 0
    
    x = round( x * pow10( -thres ) ) * pow10( thres );
    
  end
  
end
