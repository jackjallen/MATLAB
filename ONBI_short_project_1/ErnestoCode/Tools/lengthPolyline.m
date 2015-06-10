function lengths = lengthPolyline( L )

  lengths = cumsum( [ 0 ; sqrt( sum( diff( L , [] , 1 ).^2 , 2 ) ) ] );

end
