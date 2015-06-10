function saving( fname , varargin )

  START_TIME = now;

  [varargin,VERBOSE] = parseargs( varargin , '-VERBOSE','-V'  , '$FORCE$',{true,false} );
  [varargin,i,WAIT]  = parseargs( varargin , '-WAIT',   '-W'  , '$DEFS$',10 );

  fname = fixname( fname );
  [p,f,e] = fileparts( fname );
  if isempty(e), e = '.mat'; end
  fname = [p,f,e];
  
  aux1 = tmpname( 'saving_**********.saving.mat' );
  aux2 = [p,f,e, '.saving' ];
  
  CLEANOUT = { onCleanup( @()safedelete(aux1) ) , onCleanup( @()safedelete(aux2) ) };
  
  str = '';
  if numel( varargin )
  str = sprintf( ',''%s''', varargin{:} );
  end
  
  str = sprintf( 'save( ''%s'' %s);' , aux1 , str );

  if VERBOSE, fprintf('(%s)... saving in   %s\n', datestr(now) , aux1 ); end
  evalin( 'caller' , str );

  if WAIT > 0
    WAIT_START = now;
    while ~isempty( dir( [p , '*.saving.mat' ] ) )  &&  etime( WAIT_START ) < WAIT
      if VERBOSE, fprintf( 'waiting ...\n' ); end
      pause(1);
    end
  end
  

  if VERBOSE, fprintf('(%s)... moving to   %s\n', datestr(now) , aux2 ); end
  movefile( aux1 , aux2  );

  if VERBOSE, fprintf('(%s)... rename to   %s\n', datestr(now) , fname ); end
  movefile( aux2 , fname );

  if VERBOSE
    [~,t] = etime( START_TIME );
    fprintf('Saving OK. Ellapsed: %s.\n', t );
  end
  
  
  function safedelete( fn )
    try, delete( fn ); end
  end
    
end
