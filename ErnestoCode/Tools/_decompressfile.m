function fn = decompressfile( fn )

  on_dir = [];
  if ispc
    on_dir = getenv('tmp');
    if isempty(on_dir)
      on_dir = getenv('temp');
    end
  elseif isunix
    on_dir = '/tmp';
  end
        
  try
    fn = gunzip( fn , on_dir );
  catch
    msg = sprintf('Cannot decompress %s', fn);
    error(msg);
  end          

  fn = fixname( fn{1} );

end
    
