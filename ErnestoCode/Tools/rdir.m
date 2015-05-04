function files = rdir( dirname , mask , depth )
%
% TODO, allow to dirname to be a cell with several dirnames
%



  if nargin < 1,  dirname = ''; end
  if nargin < 2,  mask = []; end
  if nargin < 3,  depth = Inf; end

  FS = filesep;
  if ~iscell( dirname ), dirname = { dirname }; end
  
  for d = 1:numel( dirname )
    %dirname{d} = fixname( dirname{d} );
    try
      while dirname{d}(end) == FS
        dirname{d}(end) = [];
      end
    catch
      dirname{d} = '';
    end

    files = getfiles( dirname{d} , depth );
  end
  
  if numel( dirname ) > 1
    [~,id] = unique( { files.name } ,'first' );
    files = files(id);
  end
  
  
  function nfiles = getfiles( dname  , depth )
    if isempty( dname )
      dname = '';
      nfiles = dir();
    else
      nfiles = dir( dname );
    end
    nfiles( arrayfun( @(x) x.isdir && ( strcmp(x.name,'.') || strcmp(x.name,'..') ) , nfiles ) ) = [];
    
    if ~numel( nfiles )
      nfiles = reshape( nfiles , [0 1] );
      return;
    end

    if ~isempty( dname )
      for i = 1:numel(nfiles)
        nfiles(i).name = [ dname , FS , nfiles(i).name ];
      end
    end

    directories_idx = arrayfun( @(x) x.isdir , nfiles );
    directories = nfiles( directories_idx );
    nfiles( directories_idx ) = [];

    if ~isempty( mask )
      invalids = builtin( 'cellfun','isempty',regexp( {nfiles.name} , mask , 'ONCE' ) );
      nfiles( invalids ) = [];
    end

    if depth > 0
      for idx = 1:numel( directories )
        nfiles = [ nfiles ; getfiles( directories(idx).name , depth - 1 ) ];
      end
    end
    
    if ~numel( nfiles ), nfiles = reshape( nfiles , [0 1] ); end
    
  end
    
%     i= 1;
%     while i<= numel(nfiles)
%       if nfiles(i).isdir
%         if ~strcmp(nfiles(i).name,'.') && ~strcmp(nfiles(i).name,'..')
%           getfiles( [ dirname filesep nfiles(i).name ] );
%         end
%         nfiles(i)= [];
%       elseif ~isempty(mask) && isempty( regexp( nfiles(i).name , mask ,'ONCE') )
%         nfiles(i) = [];
%       else
%         nfiles(i).name = [ dirname filesep nfiles(i).name ];
%         i= i+1;
%       end
%     end
%     if ~isempty( nfiles )
%       files = [files ; nfiles];
%     end

end
