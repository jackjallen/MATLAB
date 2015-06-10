function K = Kernel2FK( K , szX )
%{
  a = randn([100 100 100 3]);
  k = randn([20 1 20 2]);

  tic
  ak = imfilter(a,k,'same','circular','conv');
  toc

  tic
  AK = ifftn( bsxfun( @times , fftn( a ) , Kernel2FK( k , size(a) ) ) );
  toc
  max( abs(  ak(:) - AK(:) ))

  K = Kernel2FK( k , size(a) );
  tic
  AK = ifftn( bsxfun( @times ,  fftn( a ) , K ) );
  toc

  max( abs(  ak(:) - AK(:) ))
%}

%   persistent last_k
%   persistent last_szX
%   persistent last_K
  
  if nargin < 2, szX = size( K ); end
  
%   if isequal( last_k , K )  && isequal( last_szX , szX )
%     K = last_K;
%     disp('avoid fft');
%     return;
%   end
  
  
%   last_k   = K;
%   last_szX = szX;

  nd = max( ndims(K) , numel(szX) );
  szK = size(K);
  szX(end+1:nd) = 1;
  szK(end+1:nd) = 1;
  
  c(1:nd) = {':'};
  for d = find( szK > 1 )
    if szK(d) > szX(d)
%       warning('Kernel greater than signal.');
      
      if mod( szK(d) , 2)
        c{d} = floor( ( szK(d) - szX(d) )/2 ) +(1:szX(d));
      else
        c{d} =  ceil( ( szK(d) - szX(d) )/2 ) +(1:szX(d));
      end
      K = K( c{:} );

      szK(d) = szX(d);
    elseif szK(d) < szX(d)
      c{d} = szX(d);
      K(c{:}) = 0;
    end
    c{d} = ':';
    
  end

  
  K = circshift( K , -floor(szK/2) );
  K = fftn( K );

%   last_K = K;
  
end
