function InStruct_end( BLOCK )

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
  

%% a partir de aca, ya tengo entendido lo que el usuario queria hacer.

  
  fields = evalin('caller',['fieldnames(' , BLOCK.ARGS{1} , ')']);
  for ff = 1:numel(fields)
    palabra = fields{ff};
    BLOCK.commands = regexprep( BLOCK.commands , ['(\<' palabra '\>)'] , [ BLOCK.ARGS{1} '.' palabra ] );
  end
  
  name_tmp = ['InStructCode_'  num2str(floor(rand*10000), '%04d') '.m'];
  fileID = fopen(name_tmp,'w');
  for LL = 1:numel( BLOCK.commands ),   fprintf(fileID, '%s\n',   BLOCK.commands{LL} );  end
  fclose(fileID);
  
  try
      isfile(name_tmp,'fast');
      evalin('caller',name_tmp(1:end-2))
  catch LE
      delete(name_tmp)
      error(LE.message)
  end
  delete(name_tmp)


end
