function M = FixValences( M , maxITS )

  if nargin < 2, maxITS = 1000; end

  E  = Tri2Edges( M );
  Bs = BoundsNodes( M , 'fast' );
  
  val = accumarray( E(:) , 1 );
  val( Bs ) = 2*val(Bs) - 1;
  
  for it = 1:maxITS
%     it
    [v1,n1] = max( val );
    if v1 < 6, break; end

    %     if any( ismembc( Bs , n1 ) )

    ns = setdiff( unique( E( any( ismembc( E , n1 ) , 2 ) , : ) ) , n1 );

    [v2,n2] = max( val(ns) );
    if v2 > 7
      val( n1 ) = -val(n1);
      continue;
    end
    
    n2 = ns(n2);

    TS_idx = find( any( ismembc( M.tri , n1 ) , 2 )  & any( ismembc( M.tri , n2 ) , 2 ) );

    if numel( TS_idx ) ~= 2
      val( n1 ) = -val(n1);
      continue;
    end

    m1 = setdiff( unique( M.tri( TS_idx , : ) ) , [n1 n2] );
    if numel( m1 ) ~= 2
      val( n1 ) = -val(n1);
      continue;
    end
    m2 = m1(2); m1 = m1(1);

    M.tri( TS_idx , : ) = sort( [ m1 m2 n1; m1 m2 n2] , 2 );
    val( [ n1 n2 ] ) = val( [ n1 n2 ] ) - 1;
    val( [ m1 m2 ] ) = val( [ m1 m2 ] ) + 1;
    E( all( ismembc( E , unique([n1 n2]) ) , 2 ) , : ) = sort([m1 m2]);
    hist( abs( val ) , max(abs(val)) ); drawnow
  end

  
end

