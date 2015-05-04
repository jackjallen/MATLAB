function f = filenamecase( f )

  [ f , name , ext ] = fileparts( f );
  name = [ name ext ];
  f
  name
  
  if isdir( f )
    
    names = dir( f );
    names = { names.name };
    names = names( cellfun( @(n) strcmpi( n , name ) , names ) );
    if numel( names ) == 1
      f = fullfile( f , names{1} );
    end
    
  else
    
    f = fullfile( filenamecase(f) , name );
    
  end



end
