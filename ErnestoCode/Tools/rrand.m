function x = rrand( sz , min_max )
  
  if nargin < 2
    min_max = [0 1];
  end

  if numel( min_max ) ~= 2 || ~issorted( min_max )
    errror( 'min_max  is not a sorted vector of 2 elements.');
  end

  x = rand(sz);
  
  x = x * ( min_max(2) - min_max(1) ) + min_max(1);

end