function NL = numLinesFile( fid , varargin )

  if ischar(fid)
  
    fid = fopen( fid , 'r' );
    if fid < 0
      error('imposible to open file.');
    end
    CLEANUP = onCleanup( @()fclose(fid) );
  
  else
    
    try
      pos_original = ftell(fid);
      CLEANUP = onCleanup( @()fseek( fid , pos_original , 'bof' ) );
    catch
      error('no parece ser un fid');
    end
    
  end

  [varargin,i,buffersize] = parseargs( varargin ,'BufferSize','buffer' ,'$DEFS$', 1e6 );
  [varargin,i,delimiters] = parseargs( varargin ,'DELimiters'          ,'$DEFS$', uint8( sprintf('\n') ) );
  
  NL = 0;
  
  while ~feof(fid)
    data = fread( fid , buffersize , '*uint8' );
    
    NL = NL + sum( ismembc( data , delimiters ) );
  end

end
