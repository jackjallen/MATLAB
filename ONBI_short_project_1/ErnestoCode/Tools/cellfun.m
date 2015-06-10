function varargout = cellfun( F , varargin )

  try
    switch func2str(F)
      case 'numel',     F = 'prodofsize';
      case 'isempty',   F = 'isempty';
      case 'islogical', F = 'islogical';
      case 'isreal',    F = 'isreal';
      case 'length',    F = 'length';
      case 'ndims',     F = 'ndims';
    end
  end

  try,
    [varargout{1:nargout}] = builtin( 'cellfun' , F , varargin{:} );
  catch
    [varargout{1:nargout}] = builtin( 'cellfun' , F , varargin{:} ,'UniformOutput',false);
  end

end
