function e = relerr( x , y , MAX )

  e = 2*abs( ( x - y )./( x + y )  );
  
  if nargin > 2 
    if ischar( MAX )
      if strcmpi( MAX , 'max' ) || strcmpi( MAX , 'm' )
        e = max( e(:) );
      else
        error('''MAX'' string was expected.');
      end
    else
      error('A string was expected.');
    end
  end

end

