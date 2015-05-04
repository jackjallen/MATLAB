function X = DCTresize( X , newsz )

  sz = size( X );
  
  newsz( (numel(newsz)+1):numel(sz) ) = sz( numel(newsz)+1:numel(sz) );
  
  X = idctn( dctn( X ) , newsz )*sqrt( prod(newsz)/prod(sz) );




end

