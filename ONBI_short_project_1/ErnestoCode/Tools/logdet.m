function [ld,D] = logdet(A,method,bal)

  if nargin < 2 || isempty( method ), method = 'qr'; end
  if nargin < 3, bal = true; end
  method = lower(method);
  
  if strcmp( method,'xdouble' ) ||...
     strcmp( method,'xd' )      ||...
     strcmp( method,'sym' )
    bal = false;
  end
  
  if bal
    [T,P,A]=balance(A);
  end
  
  switch method
    case 'all'
      ld = [ logdet(A,'qr',0) , logdet(A,'det',0) , logdet(A,'eig',0) , logdet(A,'svd',0) , logdet(A,'chol',0) , logdet(A,'logm',0) , logdet(A,'lu',0) , logdet(A,'ldl',0) ;...
             logdet(A,'qr',1) , logdet(A,'det',1) , logdet(A,'eig',1) , logdet(A,'svd',1) , logdet(A,'chol',1) , logdet(A,'logm',1) , logdet(A,'lu',1) , logdet(A,'ldl',1) ];
    
    case 'qr'
      [Q, T] = qr(A);
      ld = sum(log(abs(diag(T))));
      
      if nargout > 1
%         D = A \ eye(size(A));
%         D = linsolve( T , eye(size(A)) , struct('UT',true) ) * Q.';
%         D = linsolve( A , eye(size(A)) ) * Q.';
%         D = T \ (Q.');
%         D = linsolve( T , Q.' , struct('UT',true) );
        Q = Q.'; D = T \ Q;
        D = D(:).';
      end
      
    case 'det'
      ld = log( det(A) );
    case 'eig'
      ld = sum(log(abs(eig(A))));
    case 'svd'
      ld = sum(log(svd(A)));
    case 'chol'
      try
        ld = 2*sum(log(diag(chol(A))));
      catch
        ld = Inf;
      end
    case 'logm'
      ld = trace( real(logm(A) ));
    case 'lu'
      ld = sum(log(abs(diag(lu(A)))));
    case 'ldl'
      [L,D,P] = ldl(A);
      ld = sum(log(diag(D)));

    case 'xdouble'
      ld = log( double( det( xdouble( A ) ) ) );
    case 'sym'
      ld = double( log( det( sym( A ) ) ) );
  
  end
  

end
