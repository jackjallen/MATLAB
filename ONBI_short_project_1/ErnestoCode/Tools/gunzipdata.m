function X = gunzipdata( Z )

  if ~isnumeric( Z )
    error('only numeric data allowed');
  end
  if ~isreal( Z )
    error('only real data allowed');
  end
  
  Z = typecast( Z(:) , 'uint8' );
  
  fname = [ tmpname '.gz'];
  fid = fopen( fname , 'w' );
  fwrite( fid , Z );
  fclose( fid );
  
  [p,d,e] = fileparts( fname );
  
  zname = gunzip( fname , p );
  
  while iscell( zname )
    zname = zname{1};
  end
  
  fid = fopen( zname , 'r' );
  X = fread( fid , Inf , '*uint8' );
  fclose( fid );
  
  delete( fname );
  delete( zname );
  
end
