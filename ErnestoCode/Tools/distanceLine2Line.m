function [d,PQ,st] = distanceLine2Line( L1 , L2 , lb , ub )
% 
% 
% [d,PQ,st] = distanceLine2Line( [x0 y0 z0;x1 y1 z1] , [x0 y0 z0;x1 y1 z1] );
% 
%{

P0 = [10 10];
P1 = [0 0];

Q0 = [2 1];
Q1 = [12 10];

[d,PQ,st] = distanceLine2Line( [P0;P1] , [Q0;Q1] , 0 , Inf );

cla;
line( [P0(1),P1(1)] , [P0(2),P1(2)] , 'color',[0 0 1] );
line( [Q0(1),Q1(1)] , [Q0(2),Q1(2)] , 'color',[0 1 0] );
line( PQ(:,1) , PQ(:,2) , 'color',[1 0 0] ,'linewidth',2);
axis equal

%}

  if nargin < 3 || isempty(lb)
    lb = [0;0];
  end
  if isscalar(lb), lb = [lb;lb]; end
  if numel(lb)~=2, error('lb???'); end
  

  if nargin < 4 || isempty(ub)
    ub = [1;1];
  end
  if isscalar(ub), ub = [ub;ub]; end
  if numel(ub)~=2, error('ub???'); end
  
  if ~all( lb(:) <= ub(:) ), error('lb > ub !!!'); end


  if size(L1,1) ~= 2, error('a [ P0 ; P1 ] matrix is expected as L1'); end
  if size(L2,1) ~= 2, error('a [ P0 ; P1 ] matrix is expected as L2'); end

  nsd = size(L1,2);
  if size(L2,2) ~= nsd, error('L1 and L2 have to be of same dim.'); end
  
  P0 = L1(1,:);  U = L1(2,:) - P0;  P0 = P0(:); U = U(:);
  Q0 = L2(1,:);  V = L2(2,:) - Q0;  Q0 = Q0(:); V = V(:);
  
  
  UU = U.'*U;
  VV = V.'*V;
  
  
  DISTANCES = [];

  %dist( mitad L1 , mitad L2 )
  DISTANCES = [ DISTANCES ; dist( 0.5 , 0.5 ) ];


  %dist( centro , L1 ) + dist( centro , L2 )
  C = mean( [ P0 , P0+U , Q0 , Q0+V ], 2 );
  s  = U.' * ( C - P0 ) / UU;
  t  = V.' * ( C - Q0 ) / VV;
  DISTANCES = [ DISTANCES ; dist( s , t ) ];
  
    
  %dist( L1 , L2 )
  W  = P0 - Q0;
  d  = U.'*W;
  e  = V.'*W;
  UV = U.'*V;
  D  = UU*VV - UV*UV;

  s = ( UV*e - VV*d ) / D;
  t = ( UU*e - UV*d ) / D;
  DISTANCES = [ DISTANCES ; dist( s , t ) ];

  
  %dist( L1 , start L2 )
  s  = U.' * ( Q0 + 0*V - P0 ) / UU;
  DISTANCES = [ DISTANCES ; dist( s , 0 ) ];

  %dist( L1 , end L2 )
  s  = U.' * ( Q0 + 1*V - P0 ) / UU;
  DISTANCES = [ DISTANCES ; dist( s , 1 ) ];

  %dist( start L1 , L2 )
  t  = V.' * ( P0 + 0*U - Q0 ) / VV;
  DISTANCES = [ DISTANCES ; dist( 0 , t ) ];

  %dist( end L1 , L2 )
  t  = V.' * ( P0 + 1*U - Q0 ) / VV;
  DISTANCES = [ DISTANCES ; dist( 1 , t ) ];

  
  DISTANCES = [ DISTANCES ; dist( 0 , 0 ) ];          %dist( start L1 , start L2 )
  DISTANCES = [ DISTANCES ; dist( 0 , 1 ) ];          %dist( start L1 , end L2   )
  DISTANCES = [ DISTANCES ; dist( 1 , 0 ) ];          %dist( end L1   , start L2 )
  DISTANCES = [ DISTANCES ; dist( 1 , 1 ) ];          %dist( end L1   , end L2   )

  DISTANCES = [ DISTANCES ; dist( lb(1) , lb(2) ) ];  %dist( lower L1 , lower L2 )
  DISTANCES = [ DISTANCES ; dist( lb(1) , ub(2) ) ];  %dist( lower L1 , upper L2 )
  DISTANCES = [ DISTANCES ; dist( ub(1) , lb(2) ) ];  %dist( upper L1 , lower L2 )
  DISTANCES = [ DISTANCES ; dist( ub(1) , ub(2) ) ];  %dist( upper L1 , upper L2 )

  %dist( L1 , lower L2 )
  if ~isinf( lb(2) )
    s  = U.' * ( Q0 + lb(2)*V - P0 ) / UU;
    DISTANCES = [ DISTANCES ; dist( s , lb(2) ) ];
  end

  %dist( L1 , upper L2 )
  if ~isinf( ub(2) )
    s  = U.' * ( Q0 + ub(2)*V - P0 ) / UU;
    DISTANCES = [ DISTANCES ; dist( s , ub(2) ) ];
  end

  %dist( lower L1 , L2 )
  if ~isinf( lb(1) )
    t  = V.' * ( P0 + lb(1)*U - Q0 ) / VV;
    DISTANCES = [ DISTANCES ; dist( lb(1) , t ) ];
  end

  %dist( upper L1 , L2 )
  if ~isinf( ub(1) )
    t  = V.' * ( P0 + ub(1)*U - Q0 ) / VV;
    DISTANCES = [ DISTANCES ; dist( ub(1) , t ) ];
  end


  DISTANCES( isinf(DISTANCES(:,1)),: ) = [];
  DISTANCES( isnan(DISTANCES(:,1)),: ) = [];
  DISTANCES( DISTANCES(:,2)<lb(1) ,: ) = [];
  DISTANCES( DISTANCES(:,2)>ub(1) ,: ) = [];
  DISTANCES( DISTANCES(:,3)<lb(2) ,: ) = [];
  DISTANCES( DISTANCES(:,3)>ub(2) ,: ) = [];


  [kk,idx] = min( DISTANCES(:,1) );
  st = DISTANCES( idx , 2:3 );


  PQ = [ ( P0 + st(1)*U ).' ; ( Q0 + st(2)*V ).' ];
  d  = sqrt( sum( diff( PQ , 1 , 1 ).^2 ) );

  
  function d = dist( s , t )
    d = [ sum( ( ( P0 + s*U ) - ( Q0 + t*V ) ).^2 ) , s , t ];
  end


end
