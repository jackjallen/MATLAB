function [cxy,d] = ClosestPointToPolyline( xy , L )

if 0

  L = [0 0;1 0;1 1;0 1;0 0];
  
  t = linspace(0.2,2*pi-0.01,120).';
  r2 = sin(5*t)+2.5;  L = [ r2.*cos(t) , r2.*sin(t) ];
  
  
  xy = randn(100,2);
  xy = ndmat( linspace(-5,5,250) , linspace(-5,5,250) );
  
  [cxy,d] = ClosestPointToPolyline( xy , L );
  
  plot( L(:,1) , L(:,2) , 'k' );
  for i=1:size(xy,1)
    hplot( [ xy(i,1) ,  cxy(i,1) ] , [ xy(i,2) ,  cxy(i,2) ] , '.-r' );
  end
  axis equal
  
end




  d = size(xy,2);
  if size(L,2) ~= d
    error('Points and polyline must be of the same dimension.');
  end
  
  nxy = size( xy ,1);
  nL  = size( L  ,1);
  
  xy = permute( xy ,[1 3 2]);
  L  = permute( L  ,[3 1 2]);
  
  D = diff( L , 1 , 2 );

  t = bsxfun( @rdivide, ...
              sum( bsxfun(@times, bsxfun(@minus, xy , L(:,1:nL-1,:) ) , D ) ,3),...
              sum( D.^2 ,3) );
  t = max(min(t,1),0);

  cxy = bsxfun(@plus, L(:,1:nL-1,:) , bsxfun(@times, t , D ) );
  d = sum( bsxfun(@minus,cxy,xy).^2 , 3);
  [d,id] = min( d , [] , 2 );
  
  t = t( sub2indv( size(t) , [ (1:nxy).' , id ] ) );
  
  L = ipermute( L ,[3 1 2]);
  D = ipermute( D ,[3 1 2]);
  
  cxy = L( id ,:) + bsxfun(@times, t , D(id,:) );
  
  
  d = sqrt(d);

end
