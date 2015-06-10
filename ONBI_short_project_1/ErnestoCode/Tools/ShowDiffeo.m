function ShowDiffeo( phi , varargin )

  if numel(varargin) && ~ischar( varargin{1} )
    every = varargin{1};
    varargin(1) = [];
  else
    every = [1 1];
  end
    
  if numel(every)==1; every=[every every]; end

  phi = double(phi);

  b= phi(:,:,1,1) + NaN;
  b( [1:every(1):end end] , : ) = 0;
  b( : , [1:every(2):end end] ) = 0;

%   figure( 'Toolbar','figure' , 'MenuBar','figure' );
  
  if ndims(phi)==4 && size(phi,4)==3
    if numel( every ) == 3
      for i = [ every(3):every(3):size(phi,3)-1 ]
        surface('xdata',phi(:,:,i,1),'ydata',phi(:,:,i,2),'zdata',phi(:,:,i,3)   ,...
                'FaceColor','none'         ,...
                'edgecolor','interp'       ,...
                'cdata',b                  ,...
                'MarkerFaceColor',[0 0 0]  ,...
                 varargin{:}               );
      end
    else
      i = round( size(phi,3)/2 );
      s =  surface('xdata',phi(:,:,i,1),'ydata',phi(:,:,i,2),'zdata',phi(:,:,i,3)   ,...
                   'FaceColor','none'         ,...
                   'edgecolor','interp'       ,...
                   'cdata',b                  ,...
                   'MarkerFaceColor',[0 0 0]  ,...
                   varargin{:}                );

      plane = uicontrol('Style','text','string',num2str(i),'Position',[102 1 20 15]);           
      control = uicontrol('Style','slider','Min',1,'Max',size(phi,3),'value',i,'Callback',@(h,e) setPLANE , 'position', [1 1 100 15] );
      try, SetAsContinuous( control ); end
      
      
      switch getfield( ver('MATLAB') , 'Version' )
        case '7.0.4'
        case '7.2'
          drawnow;
          set( gcf , 'WindowButtonMotionFcn' , ';' );
          jf= get( handle(gcf) , 'JavaFrame' );
          set( jf.fFigureClient.getClient ,'MouseWheelMovedCallback', @(j,e) MouseWheelMoved(e) );
        case '7.6'
          set( gcf ,'WindowScrollWheelFcn' , @(h,e) MouseWheelMoved(e) );
      end

    end
    view(3);
    set( gca ,'DataAspectRatio',[1 1 1]);
    
    xl = [ min(vec(phi(:,:,:,1)))  max(vec(phi(:,:,:,1))) ]; xl = (xl-mean(xl))*1.1 + mean(xl);
    yl = [ min(vec(phi(:,:,:,2)))  max(vec(phi(:,:,:,2))) ]; yl = (yl-mean(yl))*1.1 + mean(yl);
    zl = [ min(vec(phi(:,:,:,3)))  max(vec(phi(:,:,:,3))) ]; zl = (zl-mean(zl))*1.1 + mean(zl);
    
    set( gca ,'XLim', xl , 'YLim', yl , 'ZLim', zl );
    axis(gca,'vis3d');
  elseif ndims(phi)==3 && size(phi,3)==2
    
    surface('xdata',phi(:,:,1),'ydata',phi(:,:,2),'zdata',phi(:,:,2)*0   ,...
            'FaceColor','none'         ,...
            'edgecolor','interp'       ,...
            'cdata',b                  ,...
            'MarkerFaceColor',[0 0 0]  ,...
             varargin{:}               );
    set(gca,'DataAspectRatio',[1 1 1]);
  elseif ndims(phi)==3 && size(phi,3)==3
    surface('xdata',phi(:,:,1),'ydata',phi(:,:,2),'zdata',phi(:,:,3)  ,...
            'FaceColor','none'         ,...
            'edgecolor','interp'       ,...
            'cdata',b                  ,...
            'MarkerFaceColor',[0 0 0]  ,...
             varargin{:}               );
    set(gca,'DataAspectRatio',[1 1 1]);
  end

  function setPLANE
    i = round( get(control,'Value') );
    set( plane , 'String',num2str(i) );
    set( s , 'xdata',phi(:,:,i,1),'ydata',phi(:,:,i,2),'zdata',phi(:,:,i,3) );
  end
  function x=vec(x), x=x(:); end

  function SetAsContinuous( sl )
    js = handle2javaobject( sl );
    for i=1:numel(js)
      try
        jsi = handle( js{i} , 'callbackproperties');
        jsi.AdjustmentValueChangedCallback = get(sl,'Callback');
        set(sl,'Callback','');
      end
    end
  end


  function MouseWheelMoved(e)
    if hittest == s
      i = round( get(control,'Value') );
      
      try
        i = i + e.getWheelRotation * 2;
      catch
        i = i + e.VerticalScrollCount *2;
      end
      if i<1          , i=1;           end
      if i>size(phi,3), i=size(phi,3); end
      
      set(control,'Value',i);
        set( plane , 'String',num2str(i) );
        set( s , 'xdata',phi(:,:,i,1),'ydata',phi(:,:,i,2),'zdata',phi(:,:,i,3) );
    end
  end


end
