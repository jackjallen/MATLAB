function [K,M] = Matrix2Kernel( M , szX , szK )
%
% [K,R] = Matrix2Kernel( M , size_X , size_K )
% [K,R] = Matrix2Kernel( M , size_X )    
%             by default size_K == size_X (but too slow )
%
%     M*X(:) == vec( K ** X ) + R*X(:);
%
%  if R == 0 , them the transform could be written as a filter.
%
%

%{
  a = randn([20 30 5 2]);
  k = randn([4 3 3 ]);

  M =  Kernel2Matrix( k , size(a) );
  
  [K,R] = Matrix2Kernel( M , size(a) );
  
  max( abs(  K(:) - k(:) ))
  max( abs( R(:) ) )
%}

  CONV_OPTIONS = {'same', 'circular', 'conv'};

  if nargin < 3, szK = szX; end

  Xi = zeros(szX);
  N = numel(Xi);
  Xi(:) = 1:N;
  Is = Xi(:);

  
  K = zeros(szK);
  K0 = K;

  %%% Building an spiral increasing indexed array
  Ki = K;
  Ki(:) = 1:numel(Ki);
  for i = 1:ndims(Ki), ll{i} = (1:size(Ki,i))-floor(size(Ki,i)/2)-1; end
  Ki(:) = sum( ndmat( ll{:} ).^2 , 2 );
  [Ki(:),ll]= sort( Ki(:) );
  Ki(ll) = 1:numel(Ki);
  
  for i = 1:numel(Ki)
    K0( Ki==i ) = 1;
    
    Js = imfilter( Xi , K0 , CONV_OPTIONS{:} );
    Js = Js(:);

    K0( Ki==i ) = 0;
    
    f = full( sum(sum( M .* sparse(Is,Js,1,N,N) ) ) )/N;
    K( Ki==i ) = f;

    M = M - sparse(Is,Js,f,N,N);

    if max(abs(vec(M(1,:)))) < 1e-12, 
%       fprintf('sale en %d\n',i);
      break; 
    end
  end

  K = reduceKernel( K , eps(0) );
  
  function x=vec(x), x=x(:); end
end
