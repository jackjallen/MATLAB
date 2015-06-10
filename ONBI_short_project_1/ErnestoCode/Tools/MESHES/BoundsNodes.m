function ids = BoundsNodes( M , varargin )

  ids = [];
  B = vtkFeatureEdges( M , 'BoundaryEdgesOn',[],'FeatureEdgesOff',[],'NonManifoldEdgesOff',[],'ManifoldEdgesOff',[],'ColoringOff',[] );
  if isempty( B ) || ~isfield( B, 'xyz' )  || numel( B.xyz ) == 0
    return; 
  end

  ids = vtkClosestPoint( M , B.xyz );
  if numel(varargin) && strcmp( varargin{1} , 'fast')
    return;
  end
  
  
%   no_i = @(x,i) x(x~=i);
%   ids = ids( arrayfun( @(i) any(diff(reshape(sort(no_i(M.tri(any(ismembc(M.tri,i),2 ),:),i)),2,[]),1,1)) , ids ));


  ids = sort( ids );
  E = Tri2Edges( M );
  E = E( all( ismembc( E , ids ) , 2 ) , : );
  ids = chain( E );
    
  function C = chain( E )
    C = E(1,:);
    E(1,:) = [];
    while ~isempty( E )
      r = find( E(:,1) == C(end) );
      if ~isempty( r )
        C(end+1) = E(r,2); E(r,:) = [];
        continue;
      end
      r = find( E(:,2) == C(end) );
      if ~isempty( r )
        C(end+1) = E(r,1); E(r,:) = [];
        continue;
      end
      if ~isempty( E )
        C = [ C  E(1,:) ];
        E(1,:) = [];
      end
    end
  end

end
