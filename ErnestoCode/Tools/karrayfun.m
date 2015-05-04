function y = karrayfun( M , fcn , ksize )
%
%   a= rand(10,10);
%   b= karrayfun( a , @(x) median(x(:)) , [3 4] );
%   imagesc(b);
%

  if nargin<3, ksize=1; end
  ysize= 0;

  sz= size(M);
  n_dims= numel(sz);
  
  if numel(ksize)<n_dims
    ksize( n_dims     )= 1;
    ksize( ksize == 0 )= 1;
  end
  
  for d=n_dims:-1:1
    dims{d} = 1:sz(d)-ksize(d)+1;
    ksz{d}  = 1:ksize(d);
  end
  indices= ndmat( dims{:} );
  
  oones = ones(1,n_dims);
  for i=size(indices,1):-1:1
    lind = mat2cell( indices(i,:), 1 , oones );
    for d=n_dims:-1:1
      rind{d}= lind{d} - 1 + ksz{d};
    end
    if ~ysize
      try
        aux= feval( fcn , M( rind{:} ) );
        ysize= size(aux);
        for d=numel(ysize):-1:1
          ysz{d}= 1:ysize(d);
        end
      end
    end
    
    lind= [lind ysz];
    try
      y( lind{:} ) = feval( fcn , M( rind{:} ) );
    end
  end

end
