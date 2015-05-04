function K = comm(m,n)

  if nargin == 1, n = m; end
  
  %K = @(x) reshape( permute( reshape( x , [ m,n, size(x,2) ] ) , [2 1 3] ) , [m*n,size(x,2)] );

  K = sparse( 1:(m*n) , 1:(m*n) , 1 );
  
  c = reshape(1:m*n,m,n);
  c = c';
  c = c(:);
  
  K = K(c,:);


end