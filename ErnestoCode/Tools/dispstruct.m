function s = dispstruct( x , indent , varargin )

  if nargin < 2
    indent = 0;
  end

  dispsize= 0;
  
  sz = size(x);
  sz = mat2str( sz );
  sz = strrep( sz , ' ','x' );
  sz = strrep( sz , '[',['(' sprintf('%7s',class(x)) ' '] );
  sz = strrep( sz , ']','). ' );
  
  sz_str = sz;
  if ~dispsize
    sz = '';
  end

  indentation = sprintf('%*s',indent,'');
  s  = '';
  switch class(x)
    case 'char'
      if isempty(x)
        s = [ s sz '''''' ];
      else
        for i= 1:size(x,1)
          if i~=1
            s = [ s '\n' indentation ];
            sz = sz*0+32;
          end
          s = [s sz '[ ''' strrep( x(i,:) , '\' , '\\' ) ''' ]' ];
        end
      end
      
    case {'double' 'single'}
      dims= ndims( x );
      if dims == 2  && numel(x) < 100
        s = mat2str( x , 4 );
        s = strrep(s,';',' ; ');
        s = strrep(s,'[','[ ');
        s = strrep(s,']',' ]');
      else
        s = sz_str;
      end
      if length( s ) < 100
        s= [ sz s ];
      else
        s= sz;
      end
      
    case {'int8' 'int16' 'int32' 'uint8' 'uint16' 'uint32'}
      dims= ndims( x );
      if dims == 2  && numel(x) < 100
        s = mat2str( x );
        s = strrep(s,';',' ; ');
        s = strrep(s,'[','[ ');
        s = strrep(s,']',' ]');
      else
        s = sz_str;
      end
      if length( s ) < 100
        s= [ sz s ];
      else
        s= sz;
      end
      
    case 'struct'
      if numel(x)>1
        for i = 1:numel(x)
          s = [ s indentation sprintf('(%2d) ',i)  dispstruct( x(i) , indent+5 ) '\n' ];
        end
      else      
        nfields = fieldnames( x );
        max_l = max( cell2mat( cellfun(@(x) length(x),nfields,'UniformOutput',false) ) );
        for i=1:numel(nfields)
          s = [s '\n' indentation sprintf( '%*s: ', max_l , nfields{i} )];
          ss = x.(nfields{i});
          if isstruct( ss )
            s = [s sz  dispstruct(ss, indent+max_l+2) ];
          else
            s = [indentation s dispstruct(ss,indent+max_l+2,varargin{:}) ];
          end
        end
      end
      
    case 'cell'
      if numel(x) == 1
        s = [ s '   ' dispstruct( x{1} ) '\n' ];
      else
%         s = [ s '(cell) {\n' ];
        s = [ s '(cell) {  ' ];
        for i=1:numel(x)
          if i > 1
            s = [ s '  |  ' ];
          end
          s = [ s dispstruct( x{i} , indent ) ];
        end
%         s = [ s indentation '}\n' ];
        s = [ s  '  }' ];
      end
      
      
    otherwise
      s = sz;
  end
  
  if nargout<1
    try 
      if s(1:2)=='\n', s(1:2)= []; end
    end
    fprintf( [ s '\n' ] );
    clear s;
  end

end
