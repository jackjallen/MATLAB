classdef RefVar < handle
	properties
		V=[];
	end
	methods
    function R = RefVar( V )
      if nargin < 1, V = []; end
      R.V = V;
    end
    function R = setValue( R , V )
      R.V = V;
    end
    function [varargout] = subsref( R , varargin ), [varargout{1:nargout}] = subsref( R.V , varargin{:} ); end
    function [varargout] = size(    R , varargin ), [varargout{1:nargout}] = size(    R.V , varargin{:} ); end
    function [varargout] = ndims(   R , varargin ), [varargout{1:nargout}] = ndims(   R.V , varargin{:} ); end
	end
end
