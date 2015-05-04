function Z = gzipdata( X )

  if ~isnumeric( X )
    error('only numeric data allowed');
  end
  if ~isreal( X )
    error('only real data allowed');
  end

  
  X = typecast( X(:) , 'uint8' );
  
  fname = tmpname;
  fid = fopen( fname , 'w' );
  fwrite( fid , X );
  fclose( fid );
  
  [p,d,e] = fileparts( fname );
  
  zname = gzip( fname , p );
  
  while iscell( zname )
    zname = zname{1};
  end
  
  
  fid = fopen( zname , 'r' );
  Z = fread( fid , Inf , '*uint8' );
  fclose( fid );
  
  
  
  delete( fname );
  delete( zname );
  
end
