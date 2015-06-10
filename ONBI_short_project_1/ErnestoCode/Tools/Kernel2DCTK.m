function K = Kernel2DCTK( K , szX )

  if nargin < 2, szX = size( K ); end

  szX = szX + ( szX > 1 );
  szX = szX( 1:min( numel(szX) , ndims(K) ) );
  szX( size(K) == 1 ) = 1;

  K = dctn( K , szX  , '1no' );

  c = arrayfun( @(n) 1:max(1,(n-1)) , size( K ) , 'UniformOutput' , false );

  K = K( c{:} );
  
end
