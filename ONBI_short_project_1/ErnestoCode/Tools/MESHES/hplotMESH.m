function h_ = hplotMESH( varargin )
%

  ish = ishold( gca );
  
  
  hold(gca,'on');
  
  h = plotMESH( varargin{:} );
  
  if ~ish
    hold(gca,'off');
  end
  
  
  if nargout > 0
    h_ = h;
  end

end

