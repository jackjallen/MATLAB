function [ varargout ] = sdeal( varargin )
%
% sdeal( varargin )    safe deal
%

  N = max( min( nargin , nargout ) , 1 )
  varargout(1:N) = varargin(1:N);

end