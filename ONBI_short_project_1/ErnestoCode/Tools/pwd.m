function wd=pwd

  if isunix
    [result, wd] = system('pwd');
    wd = strrep( wd , char(13) , '' );
    wd = strrep( wd , char(10) , '' );
  elseif ispc
    wd = cd;
  else
    wd = cd;
  end

end
