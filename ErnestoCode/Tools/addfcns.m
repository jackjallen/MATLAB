function [ varargout ] = addfcns( Fs , varargin )

  No = max( nargout , 1 );
  
  if iscell( Fs )
    Nf = numel(Fs);
  else
    Nf = 1;
    Fs = {Fs};
  end
  
  
  for i = 1:Nf

    if iscell( Fs{i} )

      [ this_output{1:No} ] = feval( Fs{i}{2} , varargin{:} );
      for j = 1:No
        this_output{j} = Fs{i}{1} * this_output{j};
      end
      
    else

      [ this_output{1:No} ] = feval( Fs{i} , varargin{:} );

    end
    
    for j = 1:No
      if i == 1
        varargout{j} = this_output{j};
      else
        varargout{j} = this_output{j} + varargout{j};
      end
    end
    
  end
  
  
  for j = 1:No
    if ~isnumeric( varargout{j} )
      varargout{j} = getData( varargout{j} );
    end
  end
  

end
