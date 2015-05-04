function x = fopenseekreadclose( fname , seek , count , precision )

  fid = fopen( fname , 'r' );
  if fid < 0,
    x = [];
    return;
  end
  
  fseek( fid , seek , 'bof' );
  
  x = fread( fid , count , precision );
  
  fclose( fid );

end