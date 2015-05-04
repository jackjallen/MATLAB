function E = ewlr( varargin )
% 
% 
%   e = ewlr( X , Y , W , costFun , tune )
% 
%   Weighted Linear Regression
% 
%                  _ 
%     P = min   \   w_i  * c(  X_i * P   -  Y_i  )
%                  /
%                  -
%
%   where costFun   c(x) = x^2  by default
%

  
  [P,E] = wlr( varargin{:} );


end
