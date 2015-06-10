function fn = scp( fn )

  isREMOTE = regexp( fn , '^.*@.*?:.*' , 'once' );
  if isempty( isREMOTE )
    return;
  end

  [ext,ext,ext] = fileparts( fn );
  new_fn = tmpname( [ 'remote_files/remote_file_' , ext ] , -10 );


  commando = sprintf( 'scp -p -q  "%s"  "%s"' , fn , new_fn );
  [status,result] = system( commando );
  
  if status
    error('unable to copy the remote file');
  end
  
  fn = new_fn;

end
