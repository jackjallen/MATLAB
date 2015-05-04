function E = awlr( varargin )
% 
% 
%   a = awlr( X , Y , W , costFun , tune )
% 

  
  [P,E] = wlr( varargin{:} );

  E = [P; E];

end
