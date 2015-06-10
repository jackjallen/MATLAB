function vfunction_end( BLOCK )

  if nargin == 0

    BLOCK.NAME = regexprep( mfilename , '_end$' , '' );

    BLOCK = getappdata( 0 , [ 'BLOCK_'  BLOCK.NAME ] );

    switch BLOCK.PARSE_MODE
      case 'FILE'
        stack = dbstack;
        BLOCK.commands = readFile( BLOCK.FILE , BLOCK.LINE_START:stack(2).line );

      case 'DIARY'
        set( 0 , 'Diary'     , 'off'                );
        set( 0 , 'DiaryFile' , BLOCK.old_DIARYFILE  );
        set( 0 , 'Diary'     , BLOCK.old_DIARY      );
        %try, setPrompt(''); end

        fileprintf( BLOCK.diary , 'END_____BLOCK_%s_____' , BLOCK.NAME );
        fileprintf( BLOCK.diary , ' ' );

        BLOCK.commands = readFile( BLOCK.diary );
        delete( BLOCK.diary );

        start_idx = find( cellfun( @(s) strcmp( s , sprintf( 'START_____BLOCK_%s_____' , BLOCK.NAME ) ) , BLOCK.commands ) , 1 )+1;
        BLOCK.commands  = BLOCK.commands( start_idx:end );

        end_idx   = find( cellfun( @(s) strcmp( s , sprintf( 'END_____BLOCK_%s_____'   , BLOCK.NAME ) ) , BLOCK.commands ) , 1 )-1;
        BLOCK.commands  = BLOCK.commands( 1:end_idx );

        if isempty( BLOCK.commands )
          BLOCK.PARSE_MODE = 'CELL';
          BLOCK = rmfield( BLOCK , 'old_DIARYFILE' );
          BLOCK = rmfield( BLOCK , 'old_DIARY' );
          BLOCK = rmfield( BLOCK , 'diary' );

          %trying to read the current cell in editor...
          cell_focus = false;
          try,
            desktop = com.mathworks.mde.desk.MLDesktop.getInstance;
            CELL = desktop.getSelected.getComponent(0).getComponent(0).getComponent(0).getComponent(0).getComponent(1).getComponent(0).getComponent(0);
            cell_focus = CELL.hasFocus();
          end


          if cell_focus

            LINENO = CELL.getLineFromPos( CELL.getCaretPosition ) + 1;
            TEXT   = textscan( char( CELL.getText() ) , '%s' , 'Delimiter', '\n' , 'MultipleDelimsAsOne' , false );
            TEXT   = TEXT{1};


            celldivs = [ 0 ; ...
              find( ~cellfun( 'isempty' , ...
              regexp( TEXT , '^\s*%%(?:\s+|$)' , 'ONCE' ) ) ) ;...
              numel( TEXT ) + 2 ];
            celldivs = [ celldivs( find( celldivs <= LINENO , 1 ,'last' )) ; celldivs( find( celldivs > LINENO , 1 ,'first' )) ];
            TEXT = TEXT( (celldivs(1)+1):min(numel(TEXT),(celldivs(2)-1)) );
            BLOCK_commands = removeComments( TEXT );

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

            if any( BLOCK_END_idx(1:end-1)  <  BLOCK_START_idx(2:end)  )
              error( 'BLOCKPARSING:multiples_starts_ends' , 'In cell_mode, imposible to parse more than one %s' , BLOCK.NAME );
            end

            BLOCK.commands = BLOCK_commands( ( BLOCK_START_idx(1) ):( BLOCK_END_idx(end) ) );
          end

        end

    end

    rmappdata( 0 , [ 'BLOCK_'  BLOCK.NAME ] );


    BLOCK.commands = removeComments( BLOCK.commands );
    
  end

  %clc;uneval(BLOCK);
  
  if isempty( BLOCK.commands )
    warning('BLOCKPARSING:ImposibleToParse','Imposible to parse  %s ... %s_end block.' , BLOCK.NAME , BLOCK.NAME );
  else
    if0_idx = find( ~cellfun( 'isempty' , ...
                    regexp( BLOCK.commands , '^\s*if\s+(?:0|false)' ) ) , 1 , 'first' ) +1;
    BLOCK.commands = BLOCK.commands( if0_idx:end );

    end_idx = find( ~cellfun( 'isempty' , ...
                    regexp( BLOCK.commands , '^\s*end[,;\s]*$' ) ) , 1 , 'last' ) - 1;
    BLOCK.commands = BLOCK.commands( 1:end_idx );
  end
  

  
  
  
  
if 0
  BLOCK_EDT = [ 'function  '                                ,...
                 sprintf( '%s ', BLOCK.ARGS{:} )            ,...
                 sprintf( '\n\n' )                          ,...
                 sprintf( '  %s\n', BLOCK.commands{:} )     ,...
                 sprintf( '\n\n' )                          ,...
                 sprintf( 'end' )                           ,...
                 sprintf( '\n\n' )                          ];
  
  intEditor = handle(awtcreate('com.mathworks.mlwidgets.interactivecallbacks.InteractiveCallbackEditor', ...
    'Ljava.awt.Rectangle;Ljava.lang.String;Ljava.lang.String;', ...
    java.awt.Rectangle(200,400,600,400),...
    [] , BLOCK_EDT ) );
  intEditor.setVisible(true);  
end
  
  
  
  
  
  

  
  
  
  
  
  
  
  
  

  
  FCN = [ [ 'function ' BLOCK.ARGS{:} ] ; BLOCK.commands(:) ; 'end' ];
  %FCN = OneLineCode( FCN );
  
  
  palabra = @(s) strrep(strrep( s ,'\p','(?:\<\w+\>)'),'\P','(?:\s*\<\w+\>\s*)');
  
  [FCN_name,sz] = regexp( FCN{1} , palabra('function\s+(?:(?:[\[](?:\P,?)*\P?[\]]|\p)\s*[=]\s*)?(\P)') , 'tokens' , 'once', 'tokenExtents' );
  FCN_name   = FCN_name{1};
  FCN_argins = FCN{1}(sz(end)+1:end);
  FCN_args   = regexp( FCN_argins , '\w*' , 'match' );
  
  FCN_vars   = regexp( FCN(2:end-1) , '\w*' , 'match' );
  FCN_vars   = unique( setdiff( [ FCN_vars{:} ] , FCN_args ) );
  for v = 1:numel(FCN_vars)
    vn = FCN_vars{v};
    if ~isempty( regexp( vn(1) , '\d' , 'once' ) )
      FCN_vars{v} = []; continue;
    end
    
    if ~evalin('base', [ 'exist( ''' , vn , ''' , ''var'' );' ] );
      FCN_vars{v} = []; continue;
    end
    
    lineno = find( ~cellfun( 'isempty' , regexp( FCN(2:end-1) , ['\<' vn '\>'] ) ) , 1 , 'first' ) + 1;
    eqpos  = regexp( FCN{lineno} , [ '\<' vn '\>\s*=[^=]' ] , 'once' );
    if isempty( eqpos )
      continue; 
    end
    
    if ~isempty( regexp( FCN{lineno}(eqpos+1:end) , [ '\<' vn '\>' ] , 'once' ) )
      continue;
    end
    
    if ~isempty( regexp( FCN{lineno}(1:eqpos) , [ '\<' vn '\>' ] , 'once' ) )
      continue;
    end
    
    FCN_vars{v} = [];
  end
  FCN_vars = FCN_vars( ~cellfun('isempty', FCN_vars ) );
  
  

  dirname = fileparts( which( 'vfunction.m' ) );

  VFS = readFile( fullfile( dirname , 'virtual_functions.m' ) );
  already = find( ~cellfun( 'isempty' , regexpi( VFS , [ '^\s*%VFUNCTION\s*' FCN_name ]  ) ) );
  if ~isempty(already)
    nextVF = find( ~cellfun( 'isempty' , regexpi( VFS(already+1:end) , '^\s*%VFUNCTION\s*' ) ) , 1 );
    if isempty( nextVF )
      VFS = VFS(1:already-1);
    else
      VFS = VFS([1:already-1  already+nextVF:end ]);
    end
    fid = fopen( fullfile( dirname , 'virtual_functions.m' ) , 'w' );
    cellfun( @(l) fprintf(fid,'%s\n',l) , VFS );
    fclose( fid );
  end



  
  
  FCN = [ sprintf('%%VFUNCTION  %s  %%created: %s' , to30( FCN_name ) ,  datestr(now) ) ; FCN ];
  FCN = [ FCN(1:2) ; '  %variables by vfunction' ; '  %END_variables by vfunction' ; ' ' ; FCN(3:end) ; '' ];

  for v = 1:numel(FCN_vars)
    try
      vn = FCN_vars{v};
      VAR = evalin('base', vn );
      if isnumeric(VAR) && numel(VAR) < 50
        FCN = [ FCN(1:3) ; sprintf( '  %s = %s;' , to30( vn ) , uneval(VAR) ) ; FCN(4:end) ];
      else
        FCN = [ FCN(1:3) ; sprintf( '  %s = evalin(''base'',''%s'');', to30( vn ) , vn ) ; FCN(4:end) ];
      end
    end
  end
  

  if CheckCode( FCN ), error('code has errors'); end


  fid = fopen( fullfile( dirname , 'virtual_functions.m' ) , 'a' );
  fprintf( fid , '%s\n' , FCN{:} );
  fprintf( fid , '\n\n' );
  fclose(fid);
  
%   evalin( 'caller' , [ 'clear(''' , FCN_name , ''');' ] );
  evalin( 'caller', [ FCN_name  ' = @(varargin) virtual_functions(''' FCN_name ''',varargin{:});' ] );
  

  
  
  
  function str = to30( str )
    str = [ str , blanks(30-numel(str)) ];
  end
  
  
end
