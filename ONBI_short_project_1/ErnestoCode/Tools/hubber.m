function [E,dE] = hubber( x , delta )

  if nargin < 2 , delta = 0.01; end
  
  
  E = zeros(size(x));
  
  
  if isinf( delta )
    
    smallers = isfinite( x );
    E( smallers ) = x( smallers ).^2;
    
  elseif ~delta
    
    smallers = isfinite( x );
    E( smallers ) = abs( x( smallers ) );
    
  else
  
    smallers = ( abs(x) <= delta ) &  isfinite( x );
    greaters =         ~smallers   &  isfinite( x );

    E( smallers ) = x( smallers ).^2 / ( 2 * delta );
    E( greaters ) = abs( x( greaters ) ) - delta / 2;
    
  end
  

  if nargout > 1
    
    dE = zeros(size(x));

    if isinf( delta )

      dE( smallers ) = 2*x( smallers );
      
    elseif ~delta
      
      dE( smallers ) = sign( x( smallers ) );
      
    else
    
      dE( smallers ) = x( smallers )/delta;
      dE( greaters ) = sign( x( greaters ) );
      
    end
  end



end
