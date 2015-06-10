function X = num2bin( X )
  
  X = typecast( X(:) , 'uint8' );

  X = bsxfun( @bitget , X(:) , (8:-1:1) );
  
  X = X.';
  X = char( '0' + X(:).' );
  
end
