function [T,Y] = euler( fun , tspan , y , options )

  if nargin < 4 , options = odeset; end

  MaxStep = odeget( options , 'MaxStep' , ( max(tspan)-min(tspan) )/100 );
  
  if any( diff(tspan) < 0 ), error('tspan must be an increasing vector'); end
  
  
  if numel( tspan ) > 2 
    T = nan( numel(tspan) , 1          );
    Y = nan( numel( T )   , numel( y ) );
    
    T(1,1) = tspan(1);
    Y(1,:) = y(:).';
    last_Tid = 2;
  else
    T = nan( ceil( ( max(tspan)-min(tspan) )/MaxStep ) + numel(tspan) , 1          );
    Y = nan(   numel( T )                                     , numel( y ) );
    last_Tid = 1;
  end
    
  t = tspan(1);
  nextT = 2;
  while t < tspan(end)
    h = MaxStep;
    if ( t + h ) > tspan( nextT )
      h = tspan( nextT ) - t;
      nextT = nextT + 1;
    end
    
    v = feval( fun , t , y );
    
    y = y + h * v;
    t = t + h;
    
    if numel( tspan ) > 2 

      if t == tspan( nextT - 1 )
        T(last_Tid,1) = t;
        Y(last_Tid,:) = y(:).';

        last_Tid = last_Tid + 1;
      end

    else

      T(last_Tid,1) = t;
      Y(last_Tid,:) = y(:).';

      last_Tid = last_Tid + 1;

    end

  end
  
  T = T( 1:(last_Tid-1) , : );
  Y = Y( 1:(last_Tid-1) , : );

end
