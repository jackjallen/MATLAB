function x = ininterval( x , in )

  if numel(in) ~= 2
    error('expenting ininterval( y , [x0 x1] ) and return TRUE where x \in [x0 x1).');
  end
  
%   in = sort(in);
  
  x = x >= in(1) & x < in(2);

end
