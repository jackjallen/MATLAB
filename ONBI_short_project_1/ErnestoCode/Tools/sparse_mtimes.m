function C = sparse_mtimes( L , R , szL , szR , szC )
%{
L = sprand(  16 , 32 , 1 );
R = sprand(  32 , 5 , 1 );

szL = size( L ); szR = size( R ); szC = [ 16 , 5 ];

L = reshape( L , 8 , [] );
R = reshape( R , 5 , [] );

C = reshape( reshape(L,szL) * reshape(R,szR) , szC );

CC = sparse_mtimes( L , R , szL , szR , szC );

maxnorm( C - CC )
%}

  if nargin < 3, szL = size( L ); end
  if nargin < 4, szR = size( R ); end
  if nargin < 5, szC = [ szL(1) , szR(2) ]; end

  [iL,jL,vL] = find( L );
  jL = jL - 1;
  iL = iL + jL * size(L,1);
  iL = ind2subv( szL , iL );
  jL = iL(:,2);
  iL = iL(:,1);
  %vL = nonzeros( L );

  idx = true( size(R,1) , 1 );
  idx( jL ) = false;
  
  R( idx , : ) = 0;
  
  [iR,jR,vR] = find( R );
  jR = jR - 1;
  iR = iR + jR * size(R,1);
  iR = ind2subv( szR , iR );
  jR = iR(:,2);
  iR = iR(:,1);
  %vR = nonzeros( R );

  IS = []; JS = []; VS = [];
  for k = 1:numel( iR )
    ids = ( jL == iR(k) );
%     if ~any( ids ), continue; end
    if ~nnz( ids ), continue; end
    

    IS = [ IS ; iL(ids)                   ];
    JS = [ JS ; zeros(sum(ids),1) + jR(k) ];
    VS = [ VS ; vL( ids )*vR( k )         ];
  end

  IS = IS + (JS-1)*szL(1);
  IS = ind2subv( szC , IS );

  C = sparse( IS(:,1) , IS(:,2) , VS , szC(1) , szC(2) );

end
