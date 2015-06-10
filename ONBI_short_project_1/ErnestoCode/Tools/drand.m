function x = drand( varargin )
%
%
% return a double, it is a random choice in the whole 'double' domain
%

  if nargin < 1
    
    if rand > 0.5
      x = rand(1,64) > 0.5;

      x = (-1).^x(64) * realpow( 2 , x(53:63) * realpow(2,0:10).' - 1023 ) * ...
          realpow( 2 , -52 )*( 1 + x(1:52) * realpow(2,0:51).' );

    else

      x = typecast( uint32( ceil( rand(1,2)*4294967295 ) ) , 'double' );

    end

  else
    
    x = rand( varargin{:} );
    for i=1:numel(x)
      x(i) = drand();
    end
    
  end

end
