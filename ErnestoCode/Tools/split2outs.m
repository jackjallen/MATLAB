function varargout = split2outs( c )
%
% examples
%
%   [a,b]= split2outs( [1.5 2.5] + 1);
%   [a,b]= split2outs( mat2cell( magic(3),3,[1 2] ) );
%

  oldWAR= warning('query','MATLAB:mat2cell:TrailingUnityVectorArgRemoved');
  warning('off','MATLAB:mat2cell:TrailingUnityVectorArgRemoved');

  if ~iscell(c)
    c= mat2cell( c ,ones(size(c,1),1),ones(size(c,2),1),ones(size(c,3),1),ones(size(c,4),1),ones(size(c,5),1));
  end
  varargout = c;

  warning( oldWAR.state , oldWAR.identifier );

end
