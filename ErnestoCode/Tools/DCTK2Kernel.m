function K = DCTK2Kernel( K , tol )
  
  nd = ndims( K );
  sz = size(K);

  
  M = 1;
  for d = 1:nd
    if sz(d) == 1, continue; end
    
    N = sz(d);
    if size(M,1) ~= N+1
      M = 2 * cos( ( 0:N ).' * ( (0:N) * ( pi / (N) )  ) );
      M(:,[1 end]) = M(:,[1 end])/2;
      M = inv( M );
    end
    
    perm = [ d 1:d-1 d+1:nd ];  iperm( perm ) = 1:nd;
    
    K = permute( K , perm );
    K = K(:,:);
    K( N+1 , : ) = - ( M(end,1:end-1)* K ) / M( end , end );
    
    sz(d) = N + 1;
    
    K = permute( reshape( K , sz( perm ) ) , iperm );
  end
  
  K = idctn( K , [] , '1no' );
  c = arrayfun( @(n) 1:max(1,(n-1)) , size( K ) , 'UniformOutput' , false );
  K = K( c{:} );
  
  
  if nargin < 2, tol = 1e-6; end
  tol = max( abs(K(:)) )*tol;

  c(1:nd) = {':'};
  for d = 1:nd
    while size( K,d ) >= 1
      c{d} = size(K,d);
      if max(abs(vec(K(c{:})))) < tol
        K( c{:} ) = [];
      else
        break;
      end
    end
    c{d} = ':';
  end




  K = symmetricExtension( K , 'ws' );
  
end
