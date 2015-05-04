function M = Kernel2Matrix( K , szX )
%
% M = Kernel2Matrix( K , size_X )
%   Given a filter K and a 'signal' X return the Linear Operator M
%   which  M*X(:) == vec( X ** K )
%

%{ 
  Example:
  a = randn([50 40 100 3]);
  k = randn([5 3 2]);

  ak = imfilter(a,k,'same','circular','conv');
  akf = real( ifftn( bsxfun( @times , fftn(a) , Kernel2FK(k,[50 40 100]) ) ));

  M =  Kernel2Matrix( k , size(a) );
  MA = M*a(:);

  max( abs(  ak(:) - MA(:) ))
%}

  CONV_OPTIONS = {'same', 'circular', 'conv'};

  Xi = zeros(szX);
  N = numel(Xi);
  Xi(:) = 1:N;
  Is = Xi(:);

  M = sparse( [],[],[],N , N );
  K0 = K*0;
  for i = vec(find( K ~= 0 ))'
    K0(i) = 1;

    Js = imfilter( Xi , K0 , CONV_OPTIONS{:} );
    Js = Js(:);
    
    M = M + sparse( Is , Js , K(i) , N , N );

    K0(i) = 0;
  end

  function x=vec(x), x=x(:); end
end
