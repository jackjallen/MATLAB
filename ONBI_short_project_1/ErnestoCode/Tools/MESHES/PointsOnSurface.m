function PointsOnSurface( M , vn )

  if nargin == 1
    vn = 'PointsOnSurf';
  end

  XYZ = M.xyz;
  F   = M.tri;
  try
  N   = vtkComputeNormals( M , 'AutoOrientNormalsOn', [] , 'ComputePointNormalsOff', [] , 'ComputeCellNormalsOn' , [] );
  catch
    N = ComputeNormals( M );
  end
  
  hFig = figure('Renderer','openGL','rendererMode','manual');
  
  hM = patch( 'Vertices',M.xyz , 'Faces', M.tri ,'EdgeColor','none' ,'faceColor',[.7 .2 .5] );
  axis equal;
  view(3);
  lighting gouraud;
  camlight;
  
  ax = ancestortool( hM , 'axes' );

  names = uicontrol('Parent',hFig,'Style','popupmenu','String',{'NewLine'} ,'callback',@(h,e) AddLine  ,'Position',[5 10 200 14]);
  NAMES = struct( 'Name' , {} , 'Color', {} );

  drawnow;
  
  ACTIONS = OPZ_ACTIONS;
  ACTIONS(end+1)= struct('action',{{ 'PRESSBUTTON-1'  }},'FCN',  @(h,e) START );
  ACTIONS(end+1)= struct('action',{{ 'LCONTROL' 'PRESS-Z' }},'FCN',  @(h,e) UNDO );
  ACTIONS(end+1)= struct('action',{{ 'LCONTROL' 'PRESS-N' }},'FCN',  @(h,e) AddLine(1) );

  setACTIONS( hFig , ACTIONS );

  drawnow;
  
  
  function AddLine( a )
    if nargin < 1
      LineIDX = get(names,'value') - 1;
      if LineIDX > 0, return; end
    end
    
    newname = inputdlg({'New Line Name:' },'New Name.',1, { '' } );
    newname = newname{1};
    
    NAMES(end+1).Name  = newname;
    NAMES( end ).Color = rand(1,3);
    
    set( names , 'String' , [ get(names,'String') ; newname ] );
    set( names , 'value' , numel( NAMES )+1 );
  end

  
  function UNDO
    LineIDX = get(names,'value') - 1;
    if LineIDX < 1, return; end
    LineName  = NAMES(LineIDX).Name;
    
    L = findall( hFig , 'Type','line','Tag', LineName );
    if isempty(L), return; end

    xyz = [ vec( get(L,'XData') )   vec( get(L,'YData') )   vec( get(L,'ZData') ) ];
    lastNaN = find( isnan( xyz(1:end-1,1) ) , 1 , 'last' );
    
    if isempty( lastNaN ), return; end
    
    xyz = xyz(1:lastNaN,:);
    set(L,'XData',xyz(:,1),'YData',xyz(:,2),'ZData',xyz(:,3));
    saveInBase;
  end



  function START
    LineIDX = get(names,'value') - 1;
    if LineIDX < 1, return; end
    LineName  = NAMES(LineIDX).Name;
    LineColor = NAMES(LineIDX).Color;
    
    xyz = zeros(0,3);
    Ln = line( 'XData',xyz(:,1),'YData',xyz(:,2),'ZData',xyz(:,3),'Color',[0 1 0],...
      'XLimInclude','off','YLimInclude','off','ZLimInclude','off',...
      'eraseMode','none' );

    alls = 1:size(M.tri,1);
    
    STORED_STATE = setACTIONS( hFig ,'suspend');
    set( hFig , 'keyPressFcn' , ';'     );
    set( hFig , 'Pointer'     , 'circle' );

    set( hFig , 'WindowButtonUpFcn'     , @(h,e) STOP );
    set( hFig , 'WindowButtonMotionFcn' , @(h,e) AddPoint );
  
    function AddPoint
      ray = get(ax,'currentpoint');

      ZZ = ray(1,:)-ray(2,:);
      ZZ = ZZ(:);
      nz = norm(ZZ);
      if nz, ZZ = ZZ(:)/nz; end
      ZZ = null(ZZ');

      XY = XYZ*ZZ;

      % Translate vertices so that the selection point is at the origin.
      XY = bsxfun( @minus , XY , ray(1,:)*ZZ );

      valids = alls;

      sX = sum( sign( [ XY(F(valids,1),1)   XY(F(valids,2),1)   XY(F(valids,3),1) ] ) , 2 );
      valids = valids( sX ~= 3 & sX ~= -3 );

      sX = sum( sign( [ XY(F(valids,1),2)   XY(F(valids,2),2)   XY(F(valids,3),2) ] ) , 2 );
      valids = valids( sX ~= 3 & sX ~= -3 );

      s12 = sign( XY(F(valids,1),1) .* XY(F(valids,2),2) - XY(F(valids,2),1) .* XY(F(valids,1),2) );
      s23 = sign( XY(F(valids,2),1) .* XY(F(valids,3),2) - XY(F(valids,3),1) .* XY(F(valids,2),2) );
      valids = valids( s12 == s23 | ~s12 | ~s23 );

      s23 = sign( XY(F(valids,2),1) .* XY(F(valids,3),2) - XY(F(valids,3),1) .* XY(F(valids,2),2) );
      s31 = sign( XY(F(valids,3),1) .* XY(F(valids,1),2) - XY(F(valids,1),1) .* XY(F(valids,3),2) );
      valids = valids( s23 == s31 | ~s23 | ~s31 );
      
      if isempty(valids), return; end

      t = min( ( sum( XYZ( F(valids,1) , : ).* N( valids , : ) , 2 ) - N( valids , : ) * ray(1,:)' )./( N( valids , : ) * (ray(2,:)-ray(1,:))' ) );
      p = ray(1,:) + t *( ray(2,:) - ray(1,:) );
      
      xyz = [ xyz ; p ];
      
      set( Ln ,'XData',xyz(:,1),'YData',xyz(:,2),'ZData',xyz(:,3));
    end

    function STOP
      delete(Ln);
      if ~isempty(xyz)
      
        d = [ 0 ; cumsum( sqrt( sum( diff(xyz,1).^2 , 2 ) ) ) ];
        dn = unique( [ 0:.3:d(end) d(end) ] );
        xyz = Interp1D( xyz , d , dn , 'linear' );

        [kk,xyz] = vtkClosestElement( M , xyz );

        L = findall( hFig , 'Type','line','Tag',LineName );
        if isempty( L )
          L = line( 'XData',NaN,'YData',NaN,'ZData',NaN,'Linewidth',3,'Color', LineColor,...
            'Tag',LineName,'XLimInclude','off','YLimInclude','off','ZLimInclude','off' );
%           L = line( 'XData',NaN,'YData',NaN,'ZData',NaN,'LineStyle','none','marker','.',...
%             'Color', LineColor,'markerfacecolor',LineColor,...
%             'Tag',LineName,'XLimInclude','off','YLimInclude','off','ZLimInclude','off' );
        end

        set( L , 'XData' , [ vec( get(L,'XData') ) ; xyz(:,1) ; NaN ] , ...
                 'YData' , [ vec( get(L,'YData') ) ; xyz(:,2) ; NaN ] , ...
                 'ZData' , [ vec( get(L,'ZData') ) ; xyz(:,3) ; NaN ] );
        saveInBase;
      end

      setACTIONS( hFig , 'restore' , STORED_STATE );
    end
  end  

  function saveInBase
    Ls = findall( hFig , 'Type','Line' );
    for l = Ls(:)'
      tag = get(l,'Tag');
      if isempty( tag ), continue; end
      
      xyz = [ vec( get(l,'XData') )   vec( get(l,'YData') )   vec( get(l,'ZData') ) ];
      xyz( isnan(xyz(:,1)) , : ) = [];
    
      try
        assignin( 'base' , 'PointsOnSurface___' , xyz );
        evalin( 'base' , sprintf( '%s.%s = PointsOnSurface___;', vn , tag ) );
      end
    end
    try
      evalin('base','clear PointsOnSurface___;');
    end
  end

end
