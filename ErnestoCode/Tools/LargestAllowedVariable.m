function n = LargestAllowedVariable( vartype )

  if nargin < 1
    vartype = 'double';
  end

  n = memory;
  
  n = n.MaxPossibleArrayBytes - 61488;
  if n < 0
    n = 0;
    return;
  end
  
  if ischar( vartype )
    switch lower( vartype )
      case {'double'}                           , vartype = 8;
      case {'single'}                           , vartype = 4;
      case {'logical','char','int8','uint8'}    , vartype = 1;
      case {'int16','uint16'}                   , vartype = 2;
      case {'int32','uint32'}                   , vartype = 4;
      case {'int64','uint64'}                   , vartype = 8;
    end    
  end

  n = n/vartype;
  
end