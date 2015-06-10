function [A,dA_1,dA_2] = vectorAngle( X , Y )

  X = X(:);   nX = sqrt( X.' * X );
  Y = Y(:);   nY = sqrt( Y.' * Y );

  XY = ( X.' * Y ) ./ ( nX .* nY );
  A = acos( XY );
  
  if nargout > 1
    dA_1 = ( -1/sqrt( 1 - XY.^2 ) )*...
            ( Y.'  -  ( ( Y.' * X ) * X.' ./ ( nX .* nX ) ) ) ./ ( nX .* nY );
    
    %disp( relerr( dA_1 , NumericalDiff( @(z) vectorAngle( z , Y ) , X , 'i' ) )  );
  end
  if nargout > 2
    dA_2 = ( -1/sqrt( 1 - XY.^2 ) )*...
            ( X.'  -   ( ( X.' * Y ) * Y.' ) ./ ( nY .* nY )  ) ./ ( nX .* nY );
    
    %disp( relerr( dA_2 , NumericalDiff( @(z) vectorAngle( X , z ) , Y , 'i' ) )  );
  end


end