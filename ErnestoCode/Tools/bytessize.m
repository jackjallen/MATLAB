function bs = bytessize( v , units )

  if strcmpi( units , 'disk' )
    fname = [ tmpname , '.mat' ];
    save( fname , 'v' );
    bs = filesize( fname );
    delete( fname );
    return;
  end


  vname = inputname(1);
  if isempty( vname ), vname=v; end
  
  
  var = evalin('caller' , vname );
  
  bs = whos('var');
  
  bs = bs.bytes;
  
  if nargin<2
    units = 'b';
  end
  
  switch lower(units)
    case {'bit','bits'}
      bs = bs*8;
    case {'kb','kilo','kilob','kilobyte','kilobytes'}
      bs = bs/1024;
    case {'mb','mega','megab','megabyte','megabytes'}
      bs = bs/1024/1024;
  end
  
  
end

