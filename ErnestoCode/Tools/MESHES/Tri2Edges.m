function [E,L] = Tri2Edges( M )


  if isfield( M , 'vertices' ) && ~isfield( M , 'xyz' )
    M.xyz = M.vertices;
    M = rmfield( M , 'vertices' );
  end
  if isfield( M , 'faces' ) && ~isfield( M , 'tri' )
    M.tri = M.faces;
    M = rmfield( M , 'faces' );
  end

  
  M.tri = sort( M.tri , 2 );
  
  E = [ M.tri(:,[1 2]) ; M.tri(:,[2 3]) ; M.tri(:,[1 3]) ];

  ndx = 1:size(E,1);
  [ignore,ind] = sort( E( ndx , 2 ),'ascend'); ndx = ndx(ind);
  [ignore,ind] = sort( E( ndx , 1 ),'ascend'); ndx = ndx(ind);
  E = E( ndx , : );

  E = E( [ true ; ~all( ~diff( E , 1 , 1 ) , 2 ) ] , : );

  if nargout > 1
    L = sqrt( sum( ( M.xyz( E(:,1) , : ) - M.xyz( E(:,2) , : ) ).^2 , 2 ) );
  end

end
