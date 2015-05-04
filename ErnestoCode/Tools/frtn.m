function x = frtn( x , N , type )
%
% Fast Real Transform  , incluye dct, dst, idct , hartley , etc ...
%

  if nargin == 2  &&  ischar( N )
    type = N;
    N = size(x);
  else
    if nargin < 2 || isempty(N), N = size(x);  end
    if nargin < 3,               type = 'dct'; end
  end

  nd = max( ndims(x) , numel(N) ) + 1;
  sz = size( x );
   N(nd) = 1;  N(  N==0 ) = 1;
  sz(nd) = 1; sz( sz==0 ) = 1;


  cc(1:nd) = {':'};
  for d = 1:nd
    if sz(d) == N(d), continue; end
    
    if sz(d) > N(d)
      cc{d} = 1:N(d);
      x = x( cc{:} );
    else
      cc{d} = N(d);
      x( cc{:} ) = 0;
    end
      
    cc{d} = ':';
  end

  sz = size(x);
  x  = squeeze( x );
  
  switch lower(type)
    case {'dct','dctn','dct-ii','dct2e','dct2','dct-2'},         type = 'dct-ii'; 
    case {'idct','idctn','dct-iii','dct3e','dct3','dct-3'},      type = 'dct-iii'; 
    case {'dct-i','dct1e','dct1','dct-1'},                       type = 'dct-i'; 
    case {'dct-iv','dct4e','dct4','dct-4'},                      type = 'dct-iv'; 

    case {'dst-i','dst1e','dst1','dst-1'},                       type = 'dst-i'; 
    case {'dst-ii','dst2e','dst2','dst-2'},                      type = 'dst-ii'; 
    case {'dst-iii','dst3e','dst3','dst-3'},                     type = 'dst-iii'; 
    case {'dst-iv','dst4e','dst4','dst-4'},                      type = 'dst-iv'; 

    case {'h','hartley'},                                        type = 'hartley'; 

    case {'dct1no'},                                             type = 'dct1no'; 
  end

  x = frtn_mx( x , type );
  x = reshape( x , sz   );
  
end

