function L = resamplePolyline( L , every )

  lengths = cumsum( [ 0 ; sqrt( sum( diff( L , [] , 1 ).^2 , 2 ) ) ] );

  if every < 0
    every = linspace( 0 , lengths(end) , -every+1 );
  else
    every = unique( [ 0 , lengths(end) , 0:d:lengths(end) ] );
  end
  L = Interp1D( L , lengths , every );

end
