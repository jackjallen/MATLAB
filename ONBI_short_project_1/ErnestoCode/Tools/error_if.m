function r = error_if( cond , varargin )

  r = 1;
  if cond
    if nargin < 2
      varargin = { 'algun error' };
    end
    error( varargin{:} );
  end

end
