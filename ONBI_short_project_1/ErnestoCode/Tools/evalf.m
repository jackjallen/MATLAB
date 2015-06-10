function varargout = evalf( F , varargin )

  if isa( F , 'function_handle' )

    [ varargout{1:nargout} ] = builtin( 'feval' , F , varargin{:} );
    
  elseif ( isa( F , 'char' ) && numel(F)>1 && all(F(1:2) == '@(' ) || ( ~isempty( which( F ) )  &&  ~any( F == filesep ) ))
      
    [ varargout{1:nargout} ] = builtin( 'feval' , str2func( F ) , varargin{:} );

  elseif isa( F , 'char' )

      
      cwd = pwd;
      CLEANUP = onCleanup( @() cd(cwd) );
      
      [p,f,e] = fileparts( F );
      
      arroba = strfind( p , [ filesep '@' ] );
      if ~isempty( arroba )
        p = p( 1:(arroba(1)-1) );
      end
      
      cd( p );

      if ~isempty( which( f ) )

        try

          [ varargout{1:nargout} ] = builtin( 'feval' , f , varargin{:} );
          
        catch LE
          rethrow( LE )
        end

      else

        error( 'Imposible to find a function file.' );

      end
      
  else
      
      error('what are you trying to do???');
      
  end

end
