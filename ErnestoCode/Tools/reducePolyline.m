function [xy,yy] = reducePolyline( xy , N )

  if nargin < 2, error('N was expected'); end
  if size( xy , 1 ) < 3, return; end
  
  v = [NaN NaN ; diff( xy , 1 , 1 ); NaN NaN ];
  areas = abs( v(1:end-1,1).*v(2:end,2) - v(1:end-1,2).*v(2:end,1) );
  
  r = areas < 1e-14;
  xy( r , : ) = [];
  v( r , : ) = [];
  
  while size(xy,1) > N
    areas = abs( v(1:end-1,1).*v(2:end,2) - v(1:end-1,2).*v(2:end,1) );
    
    [r,r] = min( areas );
    xy( r , : ) = [];
    v( r , : ) = [];
  end
  
  
  
  
  if nargout > 1
    yy = xy(:,2);
    xy = xy(:,1);
  end

end
