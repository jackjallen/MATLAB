function vprintf( file , str , varargin )
%{

disp('llllll')
file = 'vprintf.test'; level = 1;
vprintf( {file,Inf} );
vprintf( {file,level} , 'probando vprintf' );
pause(0.1)
vprintf( {file,level} , '\rretorna el carro' );
pause(0.1)
vprintf( {file,level} , '\{retorna }rhasta retorna' );
pause(0.1)
vprintf( {file,Inf} );
pause(0.1)
vprintf( {file,level+1} , '1234567890' );
pause(0.1)
vprintf( {file,level+1} , '\5b09876' );
pause(0.1)
vprintf( {file,level+1} , '\{09876}babcde' );
pause(0.1)
vprintf( {file,Inf} );

%}


  persistent VP
  if isempty( VP )
    VP = struct( 'filename' , {'stdout'} ,'fid',{1}, 'TimerCloser',{[]} , 'PreviousText', {''},'LastLevel',0 );
  end

  if nargin < 2, str = ''; end

  level = 0;
  if iscell( file ) && numel( file ) == 2 && ischar( file{1} ) && isnumeric( file{2} )
    level = file{2};
    file  = file{1};
  elseif iscell( file ) && numel( file ) == 2 && isempty( file{1} ) && isnumeric( file{2} )
    level = file{2};
    file  = [];
  elseif isempty( file )
    level = 1;
    file  = [];
  elseif ~ischar( file )
    level = file;
    file  = [];
  else
    error('error understanding file');
  end
  

  %%understanding str
  
  %  \rtexto        ... vuelve hasta el comienzo del carro
  GoToBeginOfLine = regexp( str , '^\\r(.*)' , 'tokens' );
  if ~isempty( GoToBeginOfLine )
    str = GoToBeginOfLine{1}{1};
    GoToBeginOfLine = true;
  else
    GoToBeginOfLine = false;
  end


  GoToText = regexp( str , '^\\\{(.*)\}r(.*)' , 'tokens' );
  if ~isempty( GoToText )
    str = GoToText{1}{2};
    GoToText = GoToText{1}{1};
  else
    GoToText = '';
  end


  RemovePrevChars = regexp( str , '^\\(\d)*b(.*)' , 'tokens' );
  if ~isempty( RemovePrevChars )
    str = RemovePrevChars{1}{2};
    RemovePrevChars = str2double( RemovePrevChars{1}{1} );
  else
    RemovePrevChars = 0;
  end
    
  RemovePrevText = regexp( str , '^\\\{(.*)\}b(.*)' , 'tokens' );
  if ~isempty( RemovePrevText )
    str = RemovePrevText{1}{2};
    RemovePrevText = RemovePrevText{1}{1};
  else
    RemovePrevText = '';
  end
  %%END understanding str
  
  if ~isempty( str )
    str = sprintf( str , varargin{:} );
  end
  

  
  if level > 0

    str_original = str;
    
    stdout_id = find( cellfun( @(fn) strcmp(fn,'stdout') , { VP.filename } ) , 1 );
    if isempty( stdout_id ), error('no deberia estar aqui'); end
    
    PreviousText = VP(stdout_id).PreviousText;
    LastLevel    = VP(stdout_id).LastLevel;
    
    newPrevText = str;
    if GoToBeginOfLine
      str = [ backs( numel(PreviousText) ) , str ];
    end

    if ~isempty( GoToText )
      if strncmp( PreviousText , GoToText , numel(GoToText) )
        str = [ backs( numel(PreviousText) - numel(GoToText) ) , str ];
      else
        str = sprintf('\n%s%s', GoToText , str );
      end
      newPrevText = [ GoToText , nobacks( str ) ];
    end
    
    if RemovePrevChars
      str = [ backs( min( RemovePrevChars , numel(PreviousText) ) ) , str ];
      newPrevText = [ PreviousText( 1 : end - min( RemovePrevChars , numel(PreviousText) ) ), nobacks( str ) ];
    end
    
    if ~isempty( RemovePrevText )
      if isequal( PreviousText( max(1,end-numel(RemovePrevText)+1):end ) , RemovePrevText )
        str = [ backs( numel( RemovePrevText ) ) , str ];
        newPrevText = [ PreviousText( 1 : end - numel(RemovePrevText) ), nobacks( str ) ];
      end
    end

    
%     if isinf( level )
%       keyboard
%     end
    
    if LastLevel ~= level && ~isinf( LastLevel ) && ~isempty(PreviousText)  &&  PreviousText(end) ~= 10
      fprintf('\n');
    end
    
    if ~isempty( str )
      fprintf( str );
    end
    
    VP(stdout_id).PreviousText = newPrevText;
    VP(stdout_id).LastLevel    = level;
    
    %set( gcf , 'name',uneval( VP(stdout_id).PreviousText ) );

    str = str_original;

  end
  
  if isempty( file ), return; end
  if level == 0,      return; end
  
  
  closeDelay = 2;


  fid_id = 0;
  for i = numel( VP ):-1:1
    if isequal( VP(i).filename , file )
      fid_id = i;
      break;
    end
  end
  
  
  if ~fid_id || ~isvalid( VP(fid_id).TimerCloser )

    fid_id = numel( VP ) + 1;
    
    VP(fid_id).filename       = file;

    fid    = fopen( file , 'a' );
    VP(fid_id).fid            = fid;

    VP(fid_id).PreviousText   = '';
    VP(fid_id).LastLevel      = 0;
    VP(fid_id).TimerCloser = timer('ExecutionMode','singleShot','StartDelay',closeDelay,'TimerFcn',@(h,e) fclose_stop_delete(fid,h) );
    
  else
    
    stop( VP(fid_id).TimerCloser );

    if get( VP(fid_id).TimerCloser , 'TasksExecuted' ) ~= 0
      fid    = fopen( file , 'a' );
      VP(fid_id).fid            = fid;
    else
      fid = VP(fid_id).fid;
    end

    set( VP(fid_id).TimerCloser , 'StartDelay',closeDelay ,  'TimerFcn',@(h,e) fclose_stop_delete(fid,h) );

  end
  
  
  
  if abs( VP(fid_id).LastLevel ) == abs( level )
    PreviousText = VP(fid_id).PreviousText;
  else
    PreviousText = '';
  end
  newPrevText  = str;

  if ~isempty( str )
  
    if ~isempty( GoToText )
      newPrevText = [ GoToText , str ];
      if strncmp( PreviousText , GoToText , numel(GoToText) )
        str = [ blanks( numel(GoToText) ) , str ];
      else
        str = [ GoToText , str ];
      end
    end

    if RemovePrevChars
      str = [ PreviousText( 1 : end - min( RemovePrevChars , numel(PreviousText) ) ) , str ];
      newPrevText = str;
    end

    if ~isempty( RemovePrevText )
      if isequal( PreviousText( max(1,end-numel(RemovePrevText)+1):end ) , RemovePrevText )
        str = [ PreviousText( 1 : end - numel(RemovePrevText) ), str  ];
        newPrevText = str;
      end
    end

    str = nobacks( str );
    while numel(str) && str(end) == char(10), str(end)=[]; end

    %fprintf( fid , '%s%s\n' , blanks( 2*abs(level) ) , str );
    fprintf( fid , '%s\n' , str );
    
  end
  
  
  VP(fid_id).PreviousText   = newPrevText;
  VP(fid_id).LastLevel      = level;
  
  start( VP(fid_id).TimerCloser );
  

  function fclose_stop_delete(fid,h)
    fclose(fid);
    stop( h );
    delete( h );
  end

    
  function s = backs( n )
    s = char( repmat( 8 , 1 , n ) );
  end
  function s = nobacks( s )
    s( s == 8 ) = [];
  end

end






% function vprintf( file , str , varargin )
% % 
% % VLEVEL is the verbose level
% % clevel is the current level
% %
% %
% %   vprintf( [] , ... )             -> print
% %   vprintf( [VLEVEL,clevel] , ... ) -> print if VLEVEL >= clevel
% %   vprintf( {'filename',VLEVEL,clevel} , ... ) -> print if VLEVEL >= clevel in the screen and in the 'filename'
% % 
% 
%   persistent VP
% 
%   if isempty( VP )
%     VP = struct( 'fileName' , {} , 'fid' , {} , 'indents' , {} , 'TimerCloser' , {} , 'currCol' , {} );
%   end
%   
%   closeDelay = 2;
%   VLEVEL = Inf;
%   clevel = 0;
%   FNAME  = '';
%   fid = -1;
%   
%   if iscell( file ) && ischar( file{1} )
%     FNAME = file{1};
%     if numel( file ) == 2  && numel( file{2} ) == 2
%       VLEVEL = file{2}(1);
%       clevel = file{2}(2);
%     elseif numel( file ) == 3
%       VLEVEL = file{2};
%       clevel = file{3};
%     elseif numel( file ) == 1
%     else
%       error('not understood file');
%     end
%   elseif isnumeric( file ) || isfloat( file )
%     if numel( file ) == 2
%       VLEVEL = file(1);
%       clevel = file(2);
%     elseif numel(file) == 0
%     else
%       error('incorrect file');
%     end
%   else
%     error('bad sintax');
%   end
% 
% 
%   
%   if VLEVEL >= abs( clevel )
%     
%     returnCarriage = regexp( str , '^\\r(.*)' , 'tokens' );
%     if ~isempty( returnCarriage )
%       str = returnCarriage{1}{1};
%       returnCarriage = Inf;
%     else
%       returnCarriage = regexp( str , '^\\([+-]?\d)*b(.*)' , 'tokens' );
%       if ~isempty( returnCarriage )
%         str = returnCarriage{1}{2};
%         returnCarriage = str2double( returnCarriage{1}{1} );
%       else
%         returnCarriage = 0;
%       end
%     end
% 
% 
%     str = sprintf( str , varargin{:} );
%     
%     
%     if clevel > 0
%       if returnCarriage == 0
%         fprintf( str );
%         VP.currCol = numel( str );
%       elseif isinf( returnCarriage )
%         fprintf( char(ones(1,VP.currCol)*8) );
%         fprintf( str );
%         VP.currCol = numel( str );
%       elseif returnCarriage > 0
%         fprintf( char(ones(1,returnCarriage)*8) );
%         fprintf( str );
%         VP.currCol = VP.currCol - returnCarriage + numel( str );
%       elseif returnCarriage < 0
%         fprintf( char(ones(1,VP.currCol+returnCarriage)*8) );
%         fprintf( str );
%         VP.currCol = - returnCarriage + numel( str );
%       end
%     end
% 
%     
%     if ~isempty( FNAME )
%       thisFid = find( cellfun( @(s) isequal( s , FNAME ) , { VP.fileName } ) , 1 );
% 
%       if ~isempty( thisFid )
%         stop( VP(thisFid).TimerCloser );
%       end
% 
%       if isempty( thisFid )
%         thisFid = numel( VP ) + 1;
%         fid = fopen( FNAME , 'a' );
% 
%         VP(thisFid).fileName    = FNAME;
%         VP(thisFid).fid         = fid;
%         VP(thisFid).TimerCloser = timer('ExecutionMode','singleShot','StartDelay',closeDelay,'TimerFcn',@(varargin) fclose(fid) );
%       elseif get( VP(thisFid).TimerCloser , 'TasksExecuted' ) > 0
%         delete( VP(thisFid).TimerCloser );
%         fid = fopen( FNAME , 'a' );
% 
%         VP(thisFid).fid         = fid;
%         VP(thisFid).TimerCloser = timer('ExecutionMode','singleShot','StartDelay',closeDelay,'TimerFcn',@(varargin) fclose(fid) );
%       else
%         fid = VP(thisFid).fid;
%       end
% 
%     end
% 
%     if fid > 0
% 
%       str( str == 8 ) = [];
%       while numel(str) && str(end) == char(10), str(end)=[]; end
%       
%       if clevel >= 0
%         if isinf( returnCarriage ) || returnCarriage == 0
% 
%         elseif returnCarriage < 0
%           str = [ blanks(-returnCarriage) str ];
%         elseif returnCarriage > 0
%           try
%             str = [ blanks( VP(thisFid).currCol - returnCarriage ) str ];
%           end
%         end
% 
%         if ~isempty(str)
%           fprintf( fid , '%s\n' , str );
%           VP(thisFid).currCol = numel( str );
%         end
% 
%       elseif ~isempty(str)
% 
%         fprintf( fid , '%s\n' , str );
% 
%       end
% 
%       start( VP(thisFid).TimerCloser );
% 
%     end
%     
%   end
% 
% end
