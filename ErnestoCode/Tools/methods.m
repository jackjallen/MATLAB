function varargout = methods( this , varargin )
%  Information on class methods

  varargout{1} = {};
  if ischar( this )
    try
      varargout{1} = feval( [ this '_methods' ] , feval( this ) );
%       disp(1);
    end
  end
  
  if isempty( varargout{1} )
    [varargout{1:nargout}] = builtin('methods',this,varargin{:});
  end
  
end


%   % varargin prepare
%   if nargin == 2
%     vargin_char = ', varargin{1}';
%   else
%     vargin_char = '';
%   end
% 
%   % char to class conversion
%   if ischar(this)
%     try
%       eval(['classinstance = ' this '();']);
%       if ismethod(classinstance,'methods')
%         % call overloaded methods method
%         if nargout == 1
%           eval(['names = methods(classinstance' vargin_char ');']);
%         else
%           eval(['methods(classinstance' vargin_char ')']);
%         end
%         return;
%       end
%     catch
%     end
%   end
% 
%   % call builtin function
%   if nargout == 1
%     eval(['names = builtin(''methods'',this' vargin_char ');']);
%   else
%     eval(['builtin(''methods'',this' vargin_char ')']);
%   end
% 
% end
