function x = fconvn( x , k , varargin )
%{
  x =      zeros( 100 , 10 , 10 , 5 , 3); x(:) = 1:numel(x); x = x+100;
  k = abs( randn( 10  , 10  , 10 , 8 , 5) );

  method = 10;

  tic; s1 = imfilter( x , k , method , 'conv','same'); toc
  tic; s2 =   fconvn( x , k , method );                toc

  maxnorm( ( s1 - s2 )./( abs(s1) + abs(s2) ) )
%}


  CC = 'conv';
  [varargin,CC] = parseargs( varargin , 'conv', '$FORCE$',{'conv',CC} );
  [varargin,CC] = parseargs( varargin , 'corr', '$FORCE$',{'corr',CC} );


  BM = 'circular';
  [varargin,BM] = parseargs( varargin , 'symmetric' , '$FORCE$',{'symmetric' ,BM} );
  [varargin,BM] = parseargs( varargin , 'circular'  , '$FORCE$',{'circular'  ,BM} );
  [varargin,BM] = parseargs( varargin , 'replicate' , '$FORCE$',{'replicate' ,BM} );

  
  SZ = 'same';
  [varargin,SZ] = parseargs( varargin , 'same'  , '$FORCE$',{'same'  ,SZ} );
  [varargin,SZ] = parseargs( varargin , 'full'  , '$FORCE$',{'full'  ,SZ} );
  [varargin,SZ] = parseargs( varargin , 'valid' , '$FORCE$',{'valid' ,SZ} );

  
  value = Inf;
  if ~isempty( varargin )
    if numel( varargin ) == 1 && isnumeric( varargin{1} )
      BM = 'value';
      value = varargin{1};
    else
      error('quedan cosas en varargin');
    end
  end
  
  
  nd  = max( ndims(x) , ndims(k) );
  szX = size(x); szX(end+1:nd) = 1;
  szK = size(k); szK(end+1:nd) = 1;

  
  
  if strcmp(CC,'corr')
    for d = find( szK > 1 ), k = flipdim(k,d); end
  end
  
%   k = reduceKernel( k );

  
  idxs = arrayfun( @(sz) 1:sz , szX , 'UniformOutput', false );

  padsize = zeros(2,nd);
  for d = 1:nd
    if szK(d) == 1, continue; end
    if strcmp( SZ , 'same' ) && strcmp( BM , 'circular' ) && szK(d) <= szX(d), continue; end
    if strcmp( SZ , 'same' ) && szX(d) == 1 && value == 0, continue; end
    
    switch SZ
      case 'same'
        n1 = floor( ( szK(d)-1 )/2 );
        n2 =  ceil( ( szK(d)-1 )/2 );
        idxs{d} = ( 1:szX(d) ) + n1;
      case 'valid'
        n1 = 0;
        n2 = 0;
        idxs{d} = ceil( szK(d)/2 )-1 + ( 1:max( szX(d) - szK(d) + 1 , 0 ) );
      case 'full'
        if strcmp( value , 0 )
          n1 = floor( szK(d)/2 );
          n2 =  ceil( szK(d)/2 )-1;
          idxs{d} = 1:( szX(d) + szK(d) - 1 );
        else
          n1 = szK(d)-1;
          n2 = szK(d)-1;
          idxs{d} = ceil( szK(d)/2 ) - 1 + ( 1:( szX(d) + szK(d) - 1 ) );
        end
    end    
    
    padsize(:,d) = [n1 ; n2];

  end
  
  if any( padsize(:) )
    x = padding( x , padsize , BM , value );
  end

  szX = size(x); szX(end+1:nd) = 1;
  szX( size(k) == 1 ) = 1;

%   x = real( ifftn( bsxfun( @times , fftn( x ) , Kernel2FK( k , szX ) )));

  k = Kernel2FK( k , szX );
  x = fftn( x );

  x = bsxfun( @times , x , k );
    
  x = ifftn( x );
  x = real( x );

  x = x( idxs{:} );


end
