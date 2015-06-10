function vfunction( varargin )

  dirname = fileparts( which( 'vfunction.m' ) );
  if ~isfile( [ dirname filesep 'virtual_functions.m' ] )
    fid = fopen( [ dirname filesep 'virtual_functions.m' ] , 'w' );
    fprintf(fid,'function varargout = virtual_functions( fn , varargin ), [ varargout{1:nargout} ] = feval( fn , varargin{:} ); end\n\n\n');
    fclose(fid);
  end  
  
  if nargin == 0
    error('vfunction:no_args','invalid call to vfunction');
  end
  
  if       strcmp( varargin{1} , '$list$' ) ||  strcmp( varargin{1} , '$list' )

    
    VFS = readFile( fullfile( dirname , 'virtual_functions.m' ) );
    VFS = VFS( ~cellfun( 'isempty' , regexpi( VFS , '^\s*%VFUNCTION\s*' ) ) );
    cellfun( @(c) disp(c), VFS );

    
  elseif   strcmp( varargin{1} , '$type$' ) ||  strcmp( varargin{1} , '$type' )
    if numel(varargin) > 1 && ischar( varargin{2} )
      FCN_name = varargin{2};
    else
      error('a vfunction name is expected.');
    end
    
    VFS = readFile( fullfile( dirname , 'virtual_functions.m' ) );
    VFS_idxs = find( ~cellfun( 'isempty' ,  regexp( VFS , '^%VFUNCTION ' ) ) );
    VFS_idxs(end+1) = numel(VFS)+1;
    
    for i = 1:numel(VFS_idxs)-1
      reg_res = regexp( VFS( VFS_idxs(i) ) , [ '^%VFUNCTION\s*' FCN_name ] );
      if isempty( reg_res{1} ), continue; end
      
      cellfun( @(c) fprintf('%s\n',c), VFS( VFS_idxs(i):(VFS_idxs(i+1)-1) ) );
      break;
    end
    
    
  elseif ~isempty( regexp( varargin{end} , '.*(end[,;]?)\s*$' , 'ONCE' ) )
          %numel( varargin{end} >= 3 )  &&  isequal( varargin{end}(end-2:end) , 'end' )
    
  
    BLOCK.NAME = 'vfunction';
    
    
    palabra = @(s) strrep(strrep( s ,'\p','(?:\<\w+\>)'),'\P','(?:\s*\<\w+\>\s*)');

    FCN = [ 'function '  , varargin{:} ];
    [FCN_def , sz ]= regexp( FCN ,palabra( '(^function\s+(?:(?:[\[](?:\P,?)*\P?[\]]|\p)\s*=\s*)?\p\s*(?:[\(]\s*(?:\P,)*\P?[\)])?\s*[,; ]).*' ) , 'tokens' , 'once' , 'tokenExtents' );
    FCN_def = regexprep(  FCN_def{1} , '^function ','' );

    FCN_body = FCN(sz(end)+1:end);
    [FCN_end ,sz] = regexp( FCN_body , '.*(end[,;]?)\s*' , 'tokens' ,'tokenExtents' ,'once' );
    FCN_body = strtrim( FCN_body(1:sz(1)-1) );
    

    BLOCK.ARGS = { FCN_def };
    BLOCK.commands = { [ 'vfunction ' FCN_def ] ;...
                       'if false'              ;...
                       FCN_body                ;...
                       'end'                   ;...
                       'vfunction_end'         };
    
    setappdata(0,'VIRTUAL_FUNCTION',BLOCK );
    try
      evalin( 'caller' , 'vfunction_end( getappdata(0,''VIRTUAL_FUNCTION'' ))' );
    catch LE
      rmappdata(0,'VIRTUAL_FUNCTION');
      rethrow( LE );
    end
    rmappdata(0,'VIRTUAL_FUNCTION');


  else


    BLOCK.NAME = mfilename;
    BLOCK.ARGS = varargin;

    %%PARSING_FILE  running in file
    stack = [];
    try, stack = dbstack('-completenames'); end
    if numel( stack ) > 1
      BLOCK.PARSE_MODE = 'FILE';
      BLOCK.FILE       = stack(2).file;
      BLOCK.LINE_START = stack(2).line;

      %checkeo que haya un "condor_end" mas adelante, sino, error!!!
      BLOCK_commands = readFile( BLOCK.FILE , [ BLOCK.LINE_START 0 ] );  %%lee hasta el final
      BLOCK_commands = removeComments( BLOCK_commands );

      BLOCK_START_idx = find( ~cellfun( 'isempty' , regexp( BLOCK_commands , [ '^\s*\<'     BLOCK.NAME  '\>'    ] ) ) );
      BLOCK_END_idx   = find( ~cellfun( 'isempty' , regexp( BLOCK_commands , [ '^[,;\s]*\<' BLOCK.NAME '_end\>' ] ) ) );

      if numel( BLOCK_START_idx ) < numel( BLOCK_END_idx )
        error( 'BLOCKPARSING:more_starts_than_ends' , 'missing %s   ... there are more %s_ends than %ss',BLOCK.NAME,BLOCK.NAME,BLOCK.NAME );
      end
      if numel( BLOCK_START_idx ) > numel( BLOCK_END_idx )
        error( 'BLOCKPARSING:more_ends_than_starts' , 'missing %s_end  ... there are more %ss than %s_ends',BLOCK.NAME,BLOCK.NAME,BLOCK.NAME  );
      end
      if ~all( BLOCK_END_idx > BLOCK_START_idx )
        error( 'BLOCKPARSING:nested_starts_ends' , 'mix or nested %ss and %s_ends' , BLOCK.NAME,BLOCK.NAME );
      end

      setappdata( 0 , [ 'BLOCK_'  BLOCK.NAME ] , BLOCK );
      return;
    end



    %%PARSING_DIARY
    BLOCK.PARSE_MODE    = 'DIARY';
    BLOCK.old_DIARY     = get( 0 , 'Diary'     );
    BLOCK.old_DIARYFILE = get( 0 , 'DiaryFile' );

    BLOCK.diary = tmpname;

    fileprintf( BLOCK.diary , 'START_____BLOCK_%s_____' , BLOCK.NAME );
    set( 0 , 'DiaryFile' , BLOCK.diary );
    set( 0 , 'Diary'     , 'on'        );

    %try, setPrompt('waiting CONDOR_END> '); end
    setappdata( 0 , [ 'BLOCK_'  BLOCK.NAME ] , BLOCK );


  end

end
