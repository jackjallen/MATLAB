function [p_,err] = invmaketransform( T , f , p , TOL )
% 
% T = maketransform( 'rxyz',[40 30 20] ,'t',[10 20 30],'s',2 );
%
%
% p = invmaketransform( T , @(p) {'t',p(1:3),'rzxz',p(4:6),'s',p(7)} );
%
% [p,err] = invmaketransform( T , @(p) {'l_affine',p(8:19),'center',p(1:3),'s',p(7),'rxyz',p(4:6),'s',p(20)} );
% 
% [p,err] = invmaketransform( T , @(p) {'center',p(1:3),'rq',p(4:6),'s',p(7)} );
%

  if nargin < 3
    p = [];
  end

  if nargin < 4
    TOL = eps(1);
  end
  
  if ~isequal( size( T ) , [4 4] ) && ~isempty( T )
    error('A 4x4 matrix expected');
  end
    
    
  if ischar( f )
    switch lower(f)
      case {'similarity','sim'}
        f = @(p) {'rxyz',p(1:3),'s',p(4),'t',p(5:7)};
        if isempty(p), p = [ 0 , 0 , 0 , det( T(1:3,1:3) )^(1/3) , T(1,4) , T(2,4) , T(3,4) ].'; end
      case {'rigid'}
        f = @(p) {'rxyz',p(1:3),'t',p(4:6)};
        if isempty(p), p = [0 , 0 , 0 , T(1,4) , T(2,4) , T(3,4) ].'; end
      case {'q'}
        f = @(p) {'q',p(1:3),'s',p(4),'t',p(5:7)};
        if isempty(p), p = [0 , 0 , 0 , det( T(1:3,1:3) )^(1/3) , T(1,4) , T(2,4) , T(3,4) ].'; end
      case {'et','t'},
        f = @(p) {'t',p(1:3)};
        if isempty(p), p = T(1:3,4); end
      case 'i',
        f = @(p) {'s',p(1)};
        if isempty(p), p = det( T(1:3,1:3) )^(1/3); end
      case 'it'
        f = @(p) {'s',p(1),'t',p(2:4)};
        if isempty(p), p = [ det( T(1:3,1:3) )^(1/3); T(1:3,4) ]; end
      case 'u',
        f = @(p) {'l_s',p(1)};
        if isempty(p), p = log( det( T(1:3,1:3) ) )/3; end
      case 'ut',
        f = @(p) {'l_s',p(1),'t',p(2:4)};
        if isempty(p), p = [ log( det( T(1:3,1:3) ) )/3 ; T(1:3,4) ]; end
      case 'f'
        f = @(p) {'s',p(1:3)};
        if isempty(p), p = svd( T(1:3,1:3) ); end
      case 'ft',
        f = @(p) {'s',p(1:3),'t',p(4:6)};
        if isempty(p), p = [ svd( T(1:3,1:3) ) ; T(1:3,4) ]; end
      case 's',
        f = @(p) {'l_s',p(1:3)};
        if isempty(p), p = log( svd( T(1:3,1:3) ) ); end
      case 'st'
        f = @(p) {'l_s',p(1:3),'t',p(4:6)};
        if isempty(p), p = [ log( svd( T(1:3,1:3) ) ); T(1:3,4) ]; end
      case 'r',
        f = @(p) {'l_xyz',p(1:3)};
        if isempty(p), p = se( real(logm( T(1:3,1:3) )) , [6;7;2] ); end
      case 'rt'
        f = @(p) {'l_xyzt',p(1:6)};
        if isempty(p), p = real( logm(T) ); p = p([7;9;2;13;14;15]); end
      case 'n'
        f = @(p) {'l_xyz',p(1:3),'s',p(4)};
        if isempty(p), p = [0;0;0; det( T(1:3,1:3) )^(1/3) ]; end
      case 'nt',
        f = @(p) {'l_xyz',p(1:3),'s',p(4),'t',p(5:7)};
        if isempty(p), p = [0;0;0; det( T(1:3,1:3) )^(1/3) ; T(1:3,4) ]; end      
      case 'm'
        f = @(p) {'l_xyzs',p(1:4)};
        if isempty(p), p = real(logm(T(1:3,1:3))); p = [ p([6;7;2]) ; trace(p)/3 ]; end
      case 'mt',
        f = @(p) {'l_xyzst',p(1:7)};
        if isempty(p), p = real(logm(T)); p = [ p([7;9;2]) ; trace(p)/3 ; p([13;14;15]) ]; end
      case 'g'
        f = @(p) {'generallinear9',p(1:9)};
        if isempty(p), p = vec( T(1:3,1:3) ); end
      case 'gt',
        f = @(p) {'generallinear',p(1:12)};
        if isempty(p), p = vec( T(1:3,:) ); end
      case 'a',
        f = @(p) {'l_affine9',p(1:9)};
        if isempty(p), p = vec( real(logm(T(1:3,1:3))) ); end
      case 'at',
        f = @(p) {'l_affine12',p(1:12)};
        if isempty(p), p = real(logm(T)); p = vec(p(1:4,:)); end
      case 'v',
        f = @(p) {'l_volumepreserving9',p(1:8)};
        if isempty(p), p = zeros(8,1); end
      case 'vt',        f = @(p) {'l_volumepreserving',p(1:11)};
        if isempty(p), p = zeros(11,1); end
      case 'p' , error('todavia no se como resolverlo!!!!');
      case 'pt', error('todavia no se como resolverlo!!!!');
    end
  end
      

  norm2 = @(x) x(:).'*x(:);

  if isempty(p)
    while 1
      try
        tt = maketransform( f , p );
        break;
      end
      p = [ p ; 1/pi^exp(1) ];
      if numel(p) > 50
        error('!!!!!!!');
      end
    end
  else
    p = p(:);
  end

  if isempty(T)
    if nargout > 1, error('empty T'); end
    p_ = p;
    return;
  end

  ids_det = unique( [ find( abs( NumericalDiff( @(p) getE_det(p) , p     , 'i' ) ) > 1e-8 ) ...
                      find( abs( NumericalDiff( @(p) getE_det(p) , p+0.1 , 'i' ) ) > 1e-8 ) ] );
  if ~isempty(ids_det)
    E = getE_det( p );
    p(ids_det) = Optimize( @(pp) getE_det( setv( p , ids_det , pp ) ) , p(ids_det)     ,'methods',{'quasinewton',50,'conjugate',50,'descendneg',1,'coordinate',1},'ls',{'quadratic','golden'} ,'noplot','verbose',0,struct('MAX_ITERATIONS',150,'MIN_ENERGY',1e-20));
%     if getE_det(p) == E
      p(ids_det) = Optimize( @(pp) getE_det( setv( p , ids_det , pp ) ) , p(ids_det)+0.1 ,'methods',{'quasinewton',50,'conjugate',50,'descendneg',1,'coordinate',1},'ls',{'quadratic','golden'} ,'noplot','verbose',0,struct('MAX_ITERATIONS',150,'MIN_ENERGY',1e-20));
%     end
  end

  ids_rot = find( abs( NumericalDiff( @(p) getE_rot(p) , p     , 'i' ) ) > 1e-8 );  ids_rot = setdiff( ids_rot , ids_det );
  if ~isempty( ids_rot )
    p(ids_rot) = Optimize( @(pp) getE_rot( setv( p , ids_rot , pp ) ) , p(ids_rot) ,'methods',{'quasinewton',50,'conjugate',50,'descendneg',1,'coordinate',1},'ls',{'quadratic','golden'} ,'noplot','verbose',0,struct('MAX_ITERATIONS',150,'MIN_ENERGY',1e-20));
  end

  ids_tra = find( abs( NumericalDiff( @(p) getE_tra(p) , p     , 'i' ) ) > 1e-8 );  ids_tra = setdiff( ids_tra , [ ids_det , ids_rot ] );
  if ~isempty( ids_tra )
    p(ids_tra) = Optimize( @(pp) getE_tra( setv( p , ids_tra , pp ) ) , p(ids_tra) ,'methods',{'quasinewton',50,'conjugate',50,'descendneg',1,'coordinate',1},'ls',{'quadratic','golden'} ,'noplot','verbose',0,struct('MAX_ITERATIONS',150,'MIN_ENERGY',1e-20));
  end
  p = round(p*1000)/1000;
  
  E = getE(p);  
  while E > TOL
    p_prev = p;
    
    p = Optimize( @(p) getE(p) , p ,'methods',{'quasinewton',200,'conjugate',200,'descendneg',1,'coordinate',1,'$BREAK$'},'ls',{'quadratic','golden'} ,'noplot','verbose',0,struct('MIN_ENERGY',1e-20));

    if isequal( p_prev , p )
      break;
    end

    E = getE(p);

    if E > TOL
      disp( getE(p) );
    end

  end
  
%   p = ExhaustiveSearch( @(p) getE(p) , p , 10 , 3 ,'verbose');


  P = f( p );
  if nargout == 0
    for i = 1:numel(P)
      if isfloat( P{i} )
        for j = 1:numel( P{i} )
          P{i}(j) = eval( sprintf('%1.4e', P{i}(j) ) );
        end
      end
    end
    
    P = uneval( P );
    disp( P(2:end-1) );
  elseif nargout == 1
    p_ = p;
  elseif nargout == 2
    p_ = p;
    M = maketransform( P{:} );
    err = T - M;
  end
  
%   alpha = 1;
%   while alpha > 1e-12
%     [E,dE] = getE( p );
%     
%     disp( E );
% 
%     alpha = alpha*1.5;
%     while alpha > 1e-12
%       if  getE( p - alpha*dE(:) ) < E
%         break;
%       end
%       alpha = alpha/1.5;
%     end
%     
%     p = p - alpha*dE(:);
%   end
  

  function E = getE_det( p )

    M = maketransform( f , p );
    
    E = norm2( det( M(1:3,1:3) ) - det( T(1:3,1:3) ) );

  end


  function E = getE_rot( p )

    M = maketransform( f , p );
    
    E = norm2( M(1:3,1:3) - T(1:3,1:3) );

  end
  

  function E = getE_tra( p )

    M = maketransform( f , p );
    
    E = norm2( M(1:3,4) - T(1:3,4) );

  end


  function [E,dE] = getE( p )

    if nargout < 2
      M = maketransform( f , p );
      
      E = norm2( M - T );
    else
      [M,dM] = maketransform( f , p );
      
      E = norm2( M - T );

      dE = 2*vect( M - T )*dM;
      
      if numel(dE) ~= numel( p )
        dE = NumericalDiff( @(p) getE(p) , p , 'i' );
      end
    end
    
  end
  

end
