function figure2background( hFig )


  if nargin < 1 || isempty( hFig )
    hFig = gcf;
  end
  
  
  AXS = findall( hFig , 'type','axes' );
  for aa = 1:numel(AXS)
    saved_AXS(aa).XColor = get( AXS(aa) , 'XColor' ); set( AXS(aa) , 'XColor',get(gcf,'color') );
    saved_AXS(aa).YColor = get( AXS(aa) , 'YColor' ); set( AXS(aa) , 'YColor',get(gcf,'color') );
    saved_AXS(aa).ZColor = get( AXS(aa) , 'ZColor' ); set( AXS(aa) , 'ZColor',get(gcf,'color') );

    saved_AXS(aa).XLim = get( AXS(aa) , 'XLim' );
    saved_AXS(aa).YLim = get( AXS(aa) , 'YLim' );
    saved_AXS(aa).ZLim = get( AXS(aa) , 'ZLim' );
  end

  
  
  bmp = getframe( hFig );
  bmp = bmp.cdata;
  
  nAX = axes('Parent',hFig,'units','normalized','position',[0 0 1 1] );
  image( bmp , 'Parent' , nAX , 'XData',1:size(bmp,2),'YData',1:size(bmp,1) );
  set( nAX , 'xlim' , [0 , size(bmp,2)] + 0.5 , 'ylim' , [0 , size(bmp,1)] + 0.5 , 'visible' , 'off' );
  
  
  for aa = 1:numel(AXS)
    set( AXS(aa) , 'XColor', saved_AXS(aa).XColor );
    set( AXS(aa) , 'YColor', saved_AXS(aa).YColor );
    set( AXS(aa) , 'ZColor', saved_AXS(aa).ZColor );
    
    set( AXS(aa) , 'XLim' , saved_AXS(aa).XLim );
    set( AXS(aa) , 'YLim' , saved_AXS(aa).YLim );
    set( AXS(aa) , 'ZLim' , saved_AXS(aa).ZLim );

    set( AXS(aa) , 'color','none' );
    delete( get( AXS(aa) , 'children' ) );
  end
  
  set( hFig , 'children' , [ AXS ; nAX ] )
  
end



