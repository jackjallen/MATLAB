function [isd,dirname] = isdir( dirname )


  isd = false;

  if exist( dirname , 'dir'  ) == 7
    isd = true;
    return;
  end
  
  
  dirname = fixfilename( [ pwd  filesep   dirname ] );

  dirname = fileparts(dirname);
  while ~isempty( dirname )
    if exist( dirname , 'dir'  ) == 7
      return;
    end
    dirname = fileparts(dirname);
  end
    

end
