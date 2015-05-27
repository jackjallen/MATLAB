function [M,pID] = TidyMesh( M , th )

  if nargin < 2, th = eps(single(1)); end
  th2 = th^2;

  
  Fs = fieldnames( M );
  Fxyz = Fs( strncmp( Fs , 'xyz' , 3 ) );
  Ftri = Fs( strncmp( Fs , 'tri' , 3 ) );
  
  X = [];
  for f = Fxyz(:).'
    X = [ X , M.(f{1}) ];
  end

  %inter-nodes distance
  D = bsxfun( @minus , permute( X , [1 3 2] ) , permute( X , [3 1 2] ) );
  D = D.^2;
  D = sum( D , 3 );
  D( D > th2 ) = Inf;
  w = isfinite( D );
  D( w ) = 1:sum( w(:) );
  

  %collapse/reassign repited and closed nodes
  [~,newID] = min( D , [] , 2 );
  M.tri = newID( M.tri );


  %remove the collapsed faces
  w = all( diff( sort( M.tri , 2 ) , 1 , 2 ) , 2 );
  M.tri = M.tri(w,:);

  %remove the unused nodes and reorder the faces
  usedNodes = unique( M.tri );
  for f = Fxyz(:).'
    M.(f{1}) = M.(f{1})(usedNodes,:);
  end
  
  iusedNodes( usedNodes ) = 1:numel( usedNodes );
  M.tri = iusedNodes( M.tri );

  
  if nargin > 1
    nX = [];
    for f = Fxyz(:).'
      nX = [ nX , M.(f{1}) ];
    end
    
    D = bsxfun( @minus , permute( nX , [1 3 2] ) , permute( X , [3 1 2] ) );
    D = D.^2;
    D = sum( D , 3 );

    [~,pID] = min( D , [] , 2 );
  end
end
