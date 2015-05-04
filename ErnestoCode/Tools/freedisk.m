function f = freedisk( d )

  if nargin<1
    d = 'c:';
  end
  
  switch lower(d)
    case 'c', d='c:';
    case 'd', d='d:';
    case 'e', d='e:';
    case 'f', d='f:';
    case 'g', d='g:';
    case 'h', d='h:';
  end
  
  [s,o]=system(['dir ' d]);
  f = str2double( strrep(strrep(strrep( o( findstr(o, 'dirs  '):end ),' ',''),'dirs',''),'byteslibres','') );

end
