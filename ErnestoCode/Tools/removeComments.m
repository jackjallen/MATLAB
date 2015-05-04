function T = removeComments( T , out )

  if ischar( T ) || isnumeric( T )
    T = readFile( T );
  end
  T = T(:);

  %first, the block comments %{--%}
  idx2 = noempty( regexp( T  , '(^|\n)\s*%}\s*($|\n)' ) , 1 );
  while  ~isempty( idx2 )
    idx1 = noempty( regexp( T(1:idx2) , '(^|\n)\s*%{\s*($|\n)' ) , 1 , 'last' );
    if isempty( idx1 )
      T{idx2} = '';
    else
      T(idx1:idx2) = repmat( {''} , idx2-idx1+1 , 1 );
    end
    idx2 = noempty( regexp( T  , '(^|\n)\s*%}\s*($|\n)' ) , 1 );
  end

  idx1 = noempty( regexp( T , '(^|\n)\s*%{\s*($|\n)' ) );
  if ~isempty(idx1)
    error('removeComments:OpenBraces','The lines %s have open blocks with no closings', strrep( mat2str(idx1) , ';' , ' ') );
  end

  %now remove the comments.
%   idxs = noempty( regexp( T , '^\s*%' ) );
%   T(idxs) = repmat( {''} , numel(idxs) , 1 );
  
%   T = regexprep( T , '(^|\n)[^\S\n]*%[^\n]*' , '$1' );

  T = regexprep( T , '((^|\n)(([\]\)}\w.]''+|[^''%])+|''[^''\n]*(''''[^''\n]*)*'')*)[^\n]*' , '$1' );
  
  %remove the trailing spaces
  T = regexprep( T , '[^\S\n]+($|\n)', '$1' );

  
  
  if nargin > 1
    fid = fopen( out , 'w' );
    cellfun( @(s) fprintf( fid , '%s\n' , s ) , T );
    fclose(fid);
  end
  
  
  
  %remove the leading spaces
  T = cellfun( @(s) strtrim(s) , T , 'UniformOutput',false );
  %remove empty lines
  T = T(noempty(T));
  
  
  
  function ne = noempty( cells , varargin )
    ne = find( ~cellfun( 'isempty' , cells ) , varargin{:} );
  end
  
end
