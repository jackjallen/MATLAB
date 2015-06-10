function x = asinhm( x )

  x = logm( x + sqrtm( x*x + eye(size(x)) ) );

end
