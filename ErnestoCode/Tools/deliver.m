function varargout = deliver(varargin)

  if nargout > nargin
    error( 'deliver:TooManyOutputs' , 'too many outputs' );
  end

  varargout = varargin( 1:max( nargout , 1 ) );

end
