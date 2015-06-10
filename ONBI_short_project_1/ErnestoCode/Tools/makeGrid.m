function G = makeGrid( X , Y )

  G = [         ndmat( [ linspace( X(1) , X(end) , 1000 ) , NaN ] , Y(:).' )   ;
        fliplr( ndmat( [ linspace( Y(1) , Y(end) , 1000 ) , NaN ] , X(:).' )) ];

end