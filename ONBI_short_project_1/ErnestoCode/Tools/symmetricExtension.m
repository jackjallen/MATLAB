function x = symmetricExtension( x , type )
  
  if nargin < 2,  type = 'ws';  end

  nd = ndims( x );
  c  = arrayfun( @(n) 1:n , size(x) , 'uniformOutput' , false );
  for d = 1:nd
    switch lower( type )
      case '0ws' ,      c{d} = [                     1:size(x,d)  ,  (size(x,d)-1):-1:1  ];
      case '0hs' ,      c{d} = [                     1:size(x,d)  ,   size(x,d)   :-1:1  ];
      case 'ws'  ,      c{d} = [  size(x,d):-1:2  ,  1:size(x,d)                         ];
      case 'hs'  ,      c{d} = [  size(x,d):-1:1  ,  1:size(x,d)                         ];
      case 'wsws',      c{d} = [  size(x,d):-1:2  ,  1:size(x,d)  ,  (size(x,d)-1):-1:1  ];
      case 'hsws',      c{d} = [  size(x,d):-1:1  ,  1:size(x,d)  ,  (size(x,d)-1):-1:1  ];
      case 'wshs',      c{d} = [  size(x,d):-1:2  ,  1:size(x,d)  ,   size(x,d)   :-1:1  ];
      case 'hshs',      c{d} = [  size(x,d):-1:1  ,  1:size(x,d)  ,   size(x,d)   :-1:1  ];
      otherwise, error('invalid type');
    end
    x    = x( c{:} );
    c{d} = 1:size(x,d);
  end

end
