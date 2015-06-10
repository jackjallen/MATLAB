function x = flip( x , ds )

  if nargin < 2
    ds = 1:ndims(x);
  end
  
  c(1:ndims(x)) = {':'};
  for d = ds
    c{d} = size(x,d):-1:1;
  end
  
%   nds = setdiff( 1:ndims(x) , ds );
%   for d = nds
%     c{d} = 1:size(x,d);
%   end
  
  x = x( c{:} );
  
end
