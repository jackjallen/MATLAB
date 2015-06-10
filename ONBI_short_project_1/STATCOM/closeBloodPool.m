function M = closeBloodPool( M )

  if numel( unique( M.triRegion ) ) ~= 1
    error('a single piece must be provided');
  end
  
  Bids = find( M.xyzIsBoundary > 0 );
  Bids = Bids([1:9,26:32,18:25,10:17]);

  if 0
  plot3d( M.xyz(Bids,:) , '.' );
  text( M.xyz(Bids,1),M.xyz(Bids,2),M.xyz(Bids,3),arrayfun(@(n)sprintf('%d',n),1:numel(Bids),'un',0),'verticalalignment','bottom','horizontalalignment','left');
  patch( 'vertices', M.xyz( : , : ) , 'faces' , [ Bids , Bids([2:end 1])] );
  end
  
  C = mean( M.xyz( Bids , : ) , 1 );
  Cid = size( M.xyz , 1 ) + 1;
  
  M.xyz( Cid , : ) = C;
  M.xyzIsBoundary( Cid ) = 1;
  
  newT(:,1) = Bids;
  newT(:,2) = Bids([2:end 1]);
  newT(:,3) = Cid;
  
  M.tri = [ M.tri ; newT ];
  M.triRegion(end+1:size(M.tri,1)) = - M.triRegion(1);

  M = FixNormals( M );

end
