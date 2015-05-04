function s = strrep( s , varargin )

  FULLMODE = false;
  if islogical( varargin{end} ) || isnumeric( varargin{end} )
    FULLMODE = varargin{end};
    varargin(end) = [];
  end
  if numel( FULLMODE ) ~= 1
    error('FULLMODE must be true or false (false by default)');
  end


  while numel( varargin )
    
    if numel( varargin ) == 1
      str2 = varargin{1}; varargin(1) = [];
      str3 = '';
    else
      str2 = varargin{1}; varargin(1) = [];
      str3 = varargin{1}; varargin(1) = [];
    end
    
    if iscell( s )
      for i=1:numel(s)
        s{i} = replace_strings( s{i} , str2 , str3 );
      end
    else
      s = replace_strings( s , str2 , str3 );
    end
    
  end
  
        
  function str = replace_strings( str , str2 , str3 )
    if FULLMODE
      for it = 1:numel(str)*10
        strp = str;
        str = builtin( 'strrep' , str , str2 , str3 );
        if isequal( str , strp ), break; end
      end
      
    else
      str = builtin( 'strrep' , str , str2 , str3 );

    end
  end

end
