function [ L , N ] = connectedLabels( I , conn )
  
  if nargin < 2      ,  conn = 'minimal';                    end
  if ischar( conn )  ,  conn = conndef( ndims(I) , conn );   end
  
  L = I*0;
  for l = unique(I(:))'
    IL = bwlabeln( I==l , conn );
    L = L + ( max(L(:) ) + IL ).*( ~~IL );
  end
  
  if nargout == 2
    N(:,2) = accumarray( L(:) , 1 );
    N(:,1) = I( arrayfun( @(l) find(L(:)==l,1) , 1:size(N,1) ) );
  end

end
