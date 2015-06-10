function K = reduceKernel( K , tol )

  nd = ndims(K);

  if nargin < 2, tol = 1e-6; end
  
  if iscell( tol )

    sz = cell2mat( tol );
    
    sz( nd+1  ) = 0;
%     sz( sz==0 ) = 1;
    
    cc(1:nd) = {':'};
    for d = 1:nd
      if sz(d)==0 || isnan(sz(d)) || isinf(sz(d)), continue; end
      
      if size(K,d) ~= sz( d )
        perm = [ d  1:d-1  d+1:nd ];
        K = permute( K , perm );
        szK = size(K);

        if     size(K,1) > sz(d)

          if mod( szK(1) , 2)
            K = K( floor( ( szK(1) - sz(d) )/2 ) +(1:sz(d)) , cc{:} );
          else
            K = K( ceil( ( szK(1) - sz(d) )/2 ) +(1:sz(d)) , cc{:} );
          end
          
        elseif size(K,1) < sz(d)

          m = sz(d) - szK(1);
          
          if mod( szK(1) , 2)
            K = cat( 1 , zeros([ ceil(m/2)    szK(2:end)]) , ...
                         K , ...
                         zeros([ floor(m/2)   szK(2:end)]) );
          else
            K = cat( 1 , zeros([ floor(m/2)   szK(2:end)]) , ...
                         K , ...
                         zeros([ ceil(m/2)    szK(2:end)]) );
          end

        end

        K = ipermute( K , perm );
      end
      
    end

  elseif isscalar( tol )

    tol = max( abs(K(:)) )*tol;
    cc(1:nd) = {':'};

    for d = 1:nd
      
      if 0
        
        center = ceil( size(K,d)/2 );
        
        cc{d} = anyn( abs( K(cc{:}) ) > tol , [ 1:d-1 d+1:nd] );
        if ~any( cc{d} ), cc{d}( center ) = true; end
        
        
        c1 = find( cc{d} , 1 , 'first' );
        c2 = find( cc{d} , 1 , 'last'  );
        w = max( center - c1 , c2 - center );
        c1 = max( 1         , center-w-1 ); if ~cc{d}(c1), c1 = c1+1; end
        c2 = min( size(K,d) , center+w+1 );
        
        cc{d} = c1:c2;
        
        K = K( cc{:} );
        cc{d} = ':';
        
      else
        
        c1 = cc;  c1{d} = 1;
        c2 = cc;
        while size( K,d ) >= 1
          c2{d} = size(K,d);
          if max(abs(vec(K(c1{:})))) < tol  &&  max(abs(vec(K(c2{:})))) < tol
            K(c2{:}) = [];
            K(c1{:}) = [];
          elseif ~mod( size(K,d) , 2 )   &&     max(abs(vec(K(c1{:})))) < tol
            K(c1{:}) = [];
            break;
          else
            break;
          end
        end
        
      end

    end

  end

  
  function i = anyn( xx , ds )

    i = ~~xx;
    for dd = ds
      i = any( i , dd );
    end
    
  end

end
