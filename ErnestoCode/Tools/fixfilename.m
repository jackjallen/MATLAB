function nfname= fixfilename( varargin )

  sys= computer;

  switch sys
    case {'PCWIN','PCWIN64'}
      novalid= '*?"<>|';
    case {'DOS'}
      novalid= char([ 32:47 56:64 91:97 123:255 ]);
    case {'SOL2','HPUX','GLNX86','GLNXA64','MAC'}
      novalid= [];
  end

  rules = getappdata(0,'FIXFILENAME_RULES');
  
  n= numel( varargin );
  for i=1:n
    fn = varargin{i};
    fn = strtrim( fn );
    
    w = which( fn );
    if ~isempty( w ), fn = w; end
    
    
    SCP = regexp( fn , '(^.*@.*?):(.*)' , 'tokens' );
    if ~isempty(SCP)
      fn = SCP{1}{2};
      SCP = SCP{1}{1};
    end
    
    
%    fn = regexprep( fn , '^/extra/', '~/extra/' );
    
    
    for c = novalid
      fn( findstr(fn,c) )=[];
    end

    
    if isunix, fn = strrep( fn , '~' , getenv('HOME') ); end

    
    fn( findstr(fn,'\') )= filesep;
    fn( findstr(fn,'/') )= filesep;
    fn( findstr( fn , [filesep filesep] ) )= [];
    

    
    [pth,name,ext]= fileparts(fn);
    if isempty(pth)
      pth = pwd; 
%     else
%       pth = cd(cd(pth));
    end
    name( findstr(name,':') )= [];
     ext( findstr( ext,':') )= [];

    fn= fullfile( pth , [ name ext ] );
   
    for r = 1:size(rules,1)
      fn = strrep( fn , rules{r,1} , rules{r,2} );
    end

    if ~isempty(SCP),
      fn = [ SCP , ':' , fn ];
    end
    
    nfname{i}= fn;
  end

  if numel(nfname)==1
    nfname= nfname{1};
  end

end
