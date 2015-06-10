function InStruct( varargin )

  if nargin == 0
    error('InStruct:no_args','invalid call to InStruct');
  end
  

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
