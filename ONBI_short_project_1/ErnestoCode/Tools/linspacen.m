function y = linspacen( x1 , x2 , N )

  if nargin < 3 , N = 100; end
  
  if N < 2 || mod(N,1), error('N is expected to be an integer greater than 2'); end

%   sz1 = size(x1);
%   sz2 = size(x2);
%   
%   nd = max( numel(sz1) , numel(sz2) );
% 
%   sz1 = [ sz1 , ones(1,nd-numel(sz1)) ];
%   sz2 = [ sz2 , ones(1,nd-numel(sz2)) ];
% 
% 
%   for d = find( sz1 ~= sz2 & ( sz1 == 1 | sz2 == 1 ) )
%     if       size(x1,d) == 1  &&  size(x2,d) > 1
%       x1 = repmat( x1 , [ ones(1,d-1) , size(x2,d) , ones(1,nd) ]  );
%     elseif   size(x2,d) == 1  &&  size(x1,d) > 1
%       x2 = repmat( x2 , [ ones(1,d-1) , size(x1,d) , ones(1,nd) ]  );
%     end
%   end
% 
%   if ~isequal( size(x1) , size(x2) )
%     error('different sizes');
%   end
%   
%   sz = size(x1);
%   
%   y = zeros( [ sz N ] );


  D = bsxfun( @minus , x2 , x1 );
  
  y = bsxfun( @times , D , vec( linspace(0,1,N) , ndims(D)+1) );
  
  y = bsxfun( @plus , x1 , y );
  

end
