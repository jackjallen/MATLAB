function A = alphafuse( A , B , alpha )

  ma = isfinite( A );
  mb = isfinite( B );

  A(  ma & mb ) = alpha * A( ma & mb ) + ( 1 - alpha ) * B( ma & mb );
  A( ~ma & mb ) = B( ~ma & mb );

end
