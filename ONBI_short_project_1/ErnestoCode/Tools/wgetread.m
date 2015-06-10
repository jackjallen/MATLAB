function [str,status] = wgetread( url )

  if nargin ~= 1 || ~ischar( url )
    error('url as char expected');
  end

  fname = strrep(url,'/','_');
  fname = regexprep(fname,'[^A-Za-z0-9_]','');
  fname = [ tmpname('wgetread_') ,'___', fname(1:min(end,50)) , '.html' ];

  [a,b]=system(['wget -O ' , fname , ' ' , url ]);
  
  fid = fopen(fname,'r');
  str=fread(fid,Inf,'*char');
  fclose(fid);
  delete(fname);
  
  str = str(:).';
end
