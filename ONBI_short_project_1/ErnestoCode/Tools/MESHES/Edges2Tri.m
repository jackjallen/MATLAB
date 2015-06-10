function T = Edges2Tri( E )

  T = [];
  EE = E;
  while numel(E)
    p1 = E(1,1);
    p2 = E(1,2);
    E(1,:) = [];
    p3 = intersect( vec( EE(EE(:,1)==p1 | EE(:,2)==p1,:) ) , vec( EE(EE(:,1)==p2 | EE(:,2)==p2,:) ) )';
    p3 = setdiff( p3 ,[p1 p2] );
    
    if ~isempty(p3)
      for p = p3
        T = [ T ; p1 p2 p];
      end
%       E( E(:,1)==p1 & E(:,2)==p3 , :) = [];
%       E( E(:,1)==p2 & E(:,2)==p3 , :) = [];
    else
      T = [ T ; p1 p2 0];
      E( 1 , :) = [];
    end
  end
  
%   T(T(:,3)==0,:)= [];
  T = sort(T,2);
  T = unique( T , 'rows' );
end
