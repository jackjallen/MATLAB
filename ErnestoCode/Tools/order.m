function idx = order( varargin )

  try
    [~     ,idx]= sort( varargin{:} );
  catch
    [ignore,idx]= sort( varargin{:} );
  end


end
