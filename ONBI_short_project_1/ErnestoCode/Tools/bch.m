function [Zh,Zl] = bch( X , Y , order , commFCN )

 
  if nargin < 3 || isempty( order )
    order = 10;
  end
  if nargin < 4,
    commFCN = @(a,b) MATRIXcomm(a,b);
  end
  if isvector( X ) && size(X,2) == 1 && isvector( Y ) && size(Y,2) == 1 &&...
     isnumeric( commFCN ) && ...
     isequal( size( commFCN ) , numel(X)*[1 1 1] , numel(Y)*[1 1 1] )
    d = numel(X);
    CSTRUC = commFCN(:,:);
    commFCN = @(a,b) CSTRUCcomm(a,b);
  end
  if ~isa( commFCN , 'function_handle' )
    error('provided commFCN must be a function_handle');
  end
    
  
  try
    commFCN(X,Y);
  catch
    error('invalid commutator or invalid X and Y.');
  end
  
  HALL = []; LYNDON = []; description = [];
  load('bch_coefficients.mat');

  
  Eh = cell( find( HALL(:,1) <= order , 1 ,'last') , 1 );
  Eh{1} = X;
  Eh{2} = Y;
  Zh = Eh{1} + Eh{2};
  for k = 3:numel(Eh)
    coeff  = HALL(k,4);
    if coeff
      Zh = Zh + coeff * getHALL(k);
    end
  end

  
  El = cell( find( LYNDON(:,1) <= order , 1 ,'last') , 1 );
  El{1} = X;
  El{2} = Y;
  Zl = El{1} + El{2};
  for k = 3:numel(El)
    coeff  = LYNDON(k,4);
    if coeff
      Zl = Zl + coeff * getLYNDON(k);
    end
  end
  
  if maxnorm( Zh - Zl ) > 1e-8
    warning('too different!! HALL serie and LYNDON serie obtain different results!');
  end
  
  
  if nargout <= 1
    Zh = (Zh + Zl)/2;
  end

  
  
  function Eh_i = getHALL( i )
    if isempty( Eh{i} )
      Eh{i} = commFCN( getHALL( HALL(i,2) ) , getHALL( HALL(i,3) ) );
    end
    Eh_i = Eh{i};
  end
  function El_i = getLYNDON( i )
    if isempty( El{i} )
      El{i} = commFCN( getLYNDON( LYNDON(i,2) ) , getLYNDON( LYNDON(i,3) ) );
    end
    El_i = El{i};
  end

  function C = MATRIXcomm( A , B )
    C = A*B - B*A;
  end

  function C = CSTRUCcomm( A , B )
    C = reshape( A.' * CSTRUC , [d,d] ).'*B;
  end

end
