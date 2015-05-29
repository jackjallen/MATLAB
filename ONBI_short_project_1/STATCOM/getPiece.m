function M = getPiece( M , ps )

  F = fieldnames( M );

  w = ismember( M.triRegion , ps );
  for f = F(:).'
    if ~strncmp( f{1} , 'tri' , 3 ), continue; end
    M.(f{1}) = M.(f{1})(w,:);
  end
  
  w = unique( M.tri );
  iw = NaN( max( w ) , 1 );
  iw( w ) = 1:numel( w );
  M.tri = iw( M.tri );
  
  for f = F(:).'
    if ~strncmp( f{1} , 'xyz' , 3 ), continue; end
    M.(f{1}) = M.(f{1})(w,:);
  end

end
