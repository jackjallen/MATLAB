function S = structedit( S )

  Sorig = S;
  Sname = inputname(1);

  hFig = figure( 'Units','pixels'                                ,...
    'Toolbar','none'                                ,...
    'Renderer', 'painter'                   ,...
    'NextPlot','new'                                ,...
    'MenuBar','none'                                ,...
    'IntegerHandle','off',...
    'DoubleBuffer','off'                            ,...
    'Position',[100 120 450 400]                    ,...
    'NumberTitle','off'                             ,...
    'waitstatus', 'waiting'                  ,...
    'Name', Sname                            );


  X = 0;
  Y = 0;
  for f = fieldnames( S )'
    makeField( f );
  end

  Y = Y + 10;
  for c = get(hFig,'children')'
    set( c , 'Position' , get(c,'Position') + [100  -Y  0  0] );
  end

  scrooleablePanel( hFig );
  
  hC = get( hFig , 'Children' );
  set( hC , 'Units', 'pixels' );
  
  
  p = uicontrol('Parent', hFig , 'Style','frame','units','pixels','position',[0 0 1000 30]);
  set( p ,'foreground' , get( p ,'background' ) );
  uicontrol('Parent', hFig , 'Style','pushbutton','String','OK'    ,'Position',[10 5 60 15] , 'callback', @(h,e) callbackOK     );
  uicontrol('Parent', hFig , 'Style','pushbutton','String','Cancel','Position',[75 5 60 15] , 'callback', @(h,e) callbackCANCEL );
  
  function callbackOK
    delete( hFig );
  end
  function callbackCANCEL
    S = Sorig;
    delete( hFig );
  end
  
  set(hFig,'ResizeFcn',@(h,e) resize );
  resize;
 
  
  hVS = findall(hFig,'tag','VerticalScrool');
  set( hFig , 'WindowScrollWheelFcn' , @(hdl,ev) controlui( hVS , 'add' , - 35 * sign(ev.VerticalScrollCount) , 'o' ) );
  
  function resize
    oldUnits = get(hFig,'Units');
    set( hFig , 'Units' ,'Pixels' );
    pos = get(hFig,'Position');
    set( hFig , 'Units',oldUnits );

    pos = [1 30 pos(3)+2  pos(4)-27];
    set( hC , 'Position', pos );
  end
  

  waitfor (hFig, 'waitstatus', 'inactive');

  
  function makeField( f )

    name = Sname;
    ss = {};
    for ff = f
      ss   = [ ss   '.' ff{1} ];
      name = [ name '.' ff{1} ];
    end
    
    ss = substruct( ss{:} );
    if isstruct( subsref( S , ss ) )
      uicontrol('style','text','string',[ name '  ' ],'Position',[ X-200 Y 200 18] , 'Horizontalalignment','right','fontunits','pixels','fontsize',12);
      Y = Y - 15;
      X = X+30;
      for ff = fieldnames( subsref( S , ss ) )'
        makeField( [ f  ff{1} ] );
      end
      X = X-30;
      Y = Y - 5;
    else
      uicontrol('style','text','string',[ name ' :'],'Position',[ X-200 Y 200 18] , 'Horizontalalignment','right','fontunits','pixels','fontsize',12);
      e = uicontrol('style','edit','Position',[ X+3 Y 1000-X 20],'Horizontalalignment','left','background',[1 1 1],'fontunits','pixels','fontsize',12);
      
      switch class( subsref( S , ss ) )
        case { 'double' , 'single' , 'logical' , 'uint8' , 'int8' , 'uint16' , 'int16' , 'uint32' , 'int32' }
          str = mat2str( subsref( S , ss ) );
          str = [ '  ' str ];
          str = strrep( strrep( str , '[' , '[ ' ) , ']' , ' ]' );
          str = strrep( strrep( str , ' ' , '   ' ) , ';' , '  ;  ' );
          set( e , 'string', str );
          
          set( e , 'callback' , @(h,ev) fillData( ss , e , 'numeric' ) );

         case { 'char' }
          str = subsref( S , ss );
          set( e , 'string', str );
          
          set( e , 'callback' , @(h,ev) fillData( ss , e , 'char' ) );
       
      end
     
      Y = Y - 24;
    end
    

  end

  function fillData( ss , e , type )
    try
      switch type
        case 'numeric'
          subsasgn( S , ss , eval( get(e,'string') ) );
        case 'char'
          subsasgn( S , ss , get(e,'string') );
      end
      set( e , 'Background' , [1 1 1] );

    catch
      set( e , 'Background' , [1 0 0] );
    end   
    
  end



end