function [X,M] = NormalizeColumnsWithCDF( X , prc )

  if nargin < 2
    prc = 5:5:95;
  end
           
  M = [];           
  if size(prc,2) == 2 &&  size(prc,1) > 2
    M = prc(:,2);
    prc = prc(:,1);
  end
  prc = prc(:);
             
  if isempty( M )             
    M = zeros([numel(prc),1]);
    for j = 1:size(X,2)
      M = M + prctile( X(:,j) , prc );
    end
  
    M = M / size(X,2);
  end

  
  for j = 1:size(X,2)
    coefs = robustfit( prctile( X(:,j) , prc ) , M );
    X(:,j) = X(:,j)*coefs(2) + coefs(1);
  end

end
