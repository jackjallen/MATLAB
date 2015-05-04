function [ varargout ] = meval( x , varargin )

  N = max( min( nargin , nargout ) , 1 );
  for i = 1:N
    varargout{i} = feval( varargin{i} , x );
  end

end