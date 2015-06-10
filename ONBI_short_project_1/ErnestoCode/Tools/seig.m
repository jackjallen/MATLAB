function [v,d] = seig( M )

  if nargout < 2
    v = eig(M);

  else

    [v,d] = eig(M);

    MM = v*d*v';

    if maxnorm( MM - M ) > 1e-8
      error('no es diag!!!');
    end

    if det(v) < 0
      v(:,3) = -v(:,3);
    end

%   MM = v*d*v';
%   maxnorm( MM - M )
    
  end  
  
end
