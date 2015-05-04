function [F,E] = nonParametricFit( x , y , beta , sz , varargin )
% 
% [F,E] = nonParametricFit( x , y , beta )
% 
%     min   1/N*| F(x) - y |^2 + beta * | laplacian(F) |^2
% 
%     E(1) = 1/N*| F(x) - y |^2
%     E(2) = integral( | laplacian(F) |^2 * dx )
% 
% 
%{
x = randn(200,1);
y = 3*(x/2).^5 + 3*x + 5 + randn(size(x))*0.1;

x = randn( 200 , 1 )*10;
y = sin( x/5 )*5 + 0.1 * (x/10).^2 + 0.5 * x + randn(size(x))*1;

figure(1);
for beta = geospace( 1e12 , 1e-6 , 100 )
  F = nonParametricFit( x , y , beta , [] , 'max_time', 20 );
  plot( F ); hplot(x,y,'.r');
  pause(.1);
end
%}

if 0

x= [ ones(10,1) ; nan(2,1) ; 0 ; nan(15,1) ; 1 ; nan(5,1) ; ones(5,1) ];

 plot( x , 'mx','markersize',20,'linewidth',4 );
hplot( nonans(x,'minlaplacian' ), '-or' )

for beta = geospace( 100 , 1e-5 , 10 )
hplot( nonParametricFit( find(~isnan(x)) , nonans(x) , beta , {-10:.01:numel(x)+11} ), 'b','linewidth',1)
pause(1)
end
hplot( nonParametricFit( find(~isnan(x)) , nonans(x) , beta , {-10:.01:numel(x)+11} ), 'm','linewidth',1)


end


  if nargin < 4, sz = []; end

  if ischar( beta )
    switch lower( beta )
      case 'l'
        Lcurve = zeros(3,200);
        bb = 1;
        for beta = geospace( 1e-9 , 1e9 , size(Lcurve,2) )
          [F,E] = nonParametricFit( x , y , beta , sz , varargin{:} );
          Lcurve(1,bb) = E(1);
          Lcurve(2,bb) = E(2);
          Lcurve(3,bb) = beta;
          
          bb = bb + 1;
        end
        
        F = Lcurve;
      case 'ocv'
        F = nonParametricFit( x , y , {'ocv'}, sz , varargin{:} );
%       case 'ocvf'
%         F = nonParametricFit( x , y , {'ocvf'}, sz , varargin{:} );
      case 'docv'
        F = nonParametricFit( x , y , {'docv'}, sz , varargin{:} );
      case 'gcv'
        F = nonParametricFit( x , y , {'gcv'}, sz , varargin{:} );
      
    end  

    return;
  end


  SOLVER    = 'normal'; %'lsqr';
  NORMALIZE = false;
  init      = [];
  max_time  = Inf;
  if nargin > 4
    [varargin,SOLVER] = parseargs( varargin ,'lsqr'                 ,'$FORCE$',{'lsqr'      ,SOLVER} );
    [varargin,SOLVER] = parseargs( varargin ,'lsqrrec'              ,'$FORCE$',{'lsqrrec'   ,SOLVER} );
    [varargin,SOLVER] = parseargs( varargin ,'normal'               ,'$FORCE$',{'normal'    ,SOLVER} );
    [varargin,SOLVER] = parseargs( varargin ,'inv'                  ,'$FORCE$',{'inv'       ,SOLVER} );
    [varargin,SOLVER] = parseargs( varargin ,'pinv'                 ,'$FORCE$',{'pinv'      ,SOLVER} );
    [varargin,SOLVER] = parseargs( varargin ,'\','backslash'        ,'$FORCE$',{'backslash' ,SOLVER} );
    
    [varargin,NORMALIZE] = parseargs( varargin ,'NORMalize','$FORCE$',{true,NORMALIZE} );
    
    [varargin,i, init] = parseargs( varargin ,'init','$DEFS$',init );
    
    [varargin,i, max_time] = parseargs( varargin ,'max_time','$DEFS$', max_time );
  end

  if isvector(x), x = x(:); end
  if isvector(y), y = y(:); end
  if ~isequal( size(x,1) , size(y,1) )
    error('different sizes!!!');
  end
  
  nz = all( isfinite(x) , 2 )  &  all( isfinite(y) , 2 );
  x = x(nz,:);
  y = y(nz,:);
  
  NSD = size(x,2);
  if NSD > 3, error('NSD too large!!'); end
  N = size(x,1);
  
      
  T = size( y , 2 );
  if T > 1 && nargout > 1
    error('only one output allowed in case of multiple regressions');
  end

  
  y_mean = mean( y , 1);
  y = bsxfun( @minus , y , y_mean );
  if NORMALIZE
    y_scale = std( y , 1 , 1 );
    y = bsxfun( @rdivide , y , y_scale );
  end

  r = [ min( x , [] , 1 )  ;...
        max( x , [] , 1 )  ];
  
  switch NSD
    case 1
      border = 20;
      if isempty(sz), sz = 1000; end
      if iscell( sz )
        F   = I3D('T',1:T,'X',sz{1},'Y',0,'Z',0);
      elseif isa( sz , 'I3D' )
        F   = sz;
      else
        dx = ( r(2,1)-r(1,1) )/(sz(1)-border-1);
        gx = linspace( r(1,1)-border/2*dx , r(2,1)+border/2*dx , sz(1) );

        F   = I3D('T',1:T,'X',gx,'Y',0,'Z',0);
      end
      dx  = mean( diff( F.X ) );
      dV  = dx;

      Lx = vec( generateStencil( ( -1:1 ) * dx , 2 ) , 1 );
      Lx = imfiltermtx( zeros( size(F,1:3) ) , Lx , 'valid','conv',NaN);
      
      L = Lx;
      
    case 2
      border = 5;
      if isempty(sz), sz = 200; end
      if iscell( sz )
        F   = I3D('T',1:T,'X',sz{1},'Y',sz{2},'Z',0);
      elseif isa( sz , 'I3D' )
        F   = sz;
      else
        if numel(sz) < 2, sz(2) = sz(1); end
        
        dx = ( r(2,1)-r(1,1) )/(sz(1)-border-1);
        gx = linspace( r(1,1)-border/2*dx , r(2,1)+border/2*dx , sz(1) );

        dy = ( r(2,2)-r(1,2) )/(sz(2)-border-1);
        gy = linspace( r(1,2)-border/2*dy , r(2,2)+border/2*dy , sz(2) );

        F   = I3D('T',1:T,'X',gx,'Y',gy,'Z',0);
      end

      dx  = mean( diff( F.X ) );
      dy  = mean( diff( F.Y ) );
      dV  = dx * dy;

      
      Lx = vec( generateStencil( ( -1:1 ) * dx , 2 ) , 1 );
      Lx = imfiltermtx( zeros( size(F,1:3) ) , Lx , 'same','conv',NaN);
      
      Ly = vec( generateStencil( ( -1:1 ) * dy , 2 ) , 2 );
      Ly = imfiltermtx( zeros( size(F,1:3) ) , Ly , 'same','conv',NaN);

      L = Lx + Ly;

    case 3
      border = 5;
      if isempty(sz), sz = 100; end
      if iscell( sz )
        F   = I3D('T',1:T,'X',sz{1},'Y',sz{2},'Z',sz{3});
      elseif isa( sz , 'I3D' )
        F   = sz;
      else
        if numel(sz) < 2, sz(2) = sz(1); end
        if numel(sz) < 3, sz(3) = sz(2); end
        
        dx = ( r(2,1)-r(1,1) )/(sz(1)-border-1);
        gx = linspace( r(1,1)-border/2*dx , r(2,1)+border/2*dx , sz(1) );

        dy = ( r(2,2)-r(1,2) )/(sz(2)-border-1);
        gy = linspace( r(1,2)-border/2*dy , r(2,2)+border/2*dy , sz(2) );

        dz = ( r(2,3)-r(1,3) )/(sz(3)-border-1);
        gz = linspace( r(1,3)-border/2*dz , r(2,3)+border/2*dz , sz(3) );

        F   = I3D('T',1:T,'X',gx,'Y',gy,'Z',gz);
      end

      dx  = mean( diff( F.X ) );
      dy  = mean( diff( F.Y ) );
      dz  = mean( diff( F.Z ) );
      dV  = dx * dy * dz;

      
      Lx = vec( generateStencil( ( -1:1 ) * dx , 2 ) , 1 );
      Lx = imfiltermtx( zeros( size(F,1:3) ) , Lx , 'same','conv',NaN);
      
      Ly = vec( generateStencil( ( -1:1 ) * dy , 2 ) , 2 );
      Ly = imfiltermtx( zeros( size(F,1:3) ) , Ly , 'same','conv',NaN);

      Lz = vec( generateStencil( ( -1:1 ) * dz , 2 ) , 3 );
      Lz = imfiltermtx( zeros( size(F,1:3) ) , Lz , 'same','conv',NaN);

      L = Lx + Ly + Lz;
      
  end

  M = F.weights( resize( x , [] , 3 ) );
  
  
  if nargout > 1
    Mo = M;
    Lo = L;
    yo = y;
  end
  
  M = M / sqrt(N);
  y = y / sqrt(N);
  
  if iscell(beta) && strcmp(beta{1},'ocv')
    if diff(size(M))>=0
      f_aux = @(beta) M*((M.'*M + L.'*L*beta*dV)\M.');
    else
      f_aux = @(beta) M*((M.'*M + L.'*L*beta*dV)\eye(size(M,2)))*M.';
    end
    f = @(beta) sum(((f_aux(beta)*y - y)./(1-diag(f_aux(beta)))).^2);
    F = 10.^Optimize1D(@(b) f(10.^b) , [-7  3], [],'methods', 'golden');
    return;
%   elseif iscell(beta) && strcmp(beta{1},'ocvf')
%     f_aux = @(beta) M*(factorize((M.'*M + L.'*L*beta*dV))\M.');
%     f = @(beta) sum(((f_aux(beta)*y - y)./(1-diag(f_aux(beta)))).^2);
%     F = 10.^Optimize1D(@(b) f(10.^b) , [-7  3], [],'methods', 'golden');
%     return;
  elseif iscell(beta) && strcmp(beta{1},'gcv')
    if diff(size(M))>=0
      f_aux = @(beta) M*((M.'*M + L.'*L*beta*dV)\M.');
    else
      f_aux = @(beta) M*((M.'*M + L.'*L*beta*dV)\eye(size(M,2)))*M.';
    end
    f = @(beta) sum((f_aux(beta)*y - y).^2)/(1-trace(f_aux(beta))/N).^2;
    F = 10.^Optimize1D(@(b) f(10.^b) , [-7  3], [],'methods', 'golden');
    return;
  elseif iscell(beta) && strcmp(beta{1},'docv')
    if diff(size(M))>=0
      f_aux = @(beta) M*((M.'*M + L.'*L*beta*dV)\M.');
    else
      f_aux = @(beta) M*((M.'*M + L.'*L*beta*dV)\eye(size(M,2)))*M.';
    end
    f = @(beta) sum(((f_aux(beta)*y - y)./(1-diag(f_aux(beta)))).^2);
    F = 10.^Optimize1D(@(b) - NumericalDiff(@(x) f(x),10.^b,'i' ) , [-7  3], [],'methods', 'golden');
    return;
  end
  
  L = L * sqrt( beta * dV );
  
  switch lower( SOLVER )
    case {'normal'}
      L = ( M.'*M + L.'*L ); %cond(L)
      y = ( M.' * y );

      F.data(:) = L \ y;

    case {'inv'}
      L = ( M.'*M + L.'*L );
      y = ( M.' * y );
      
      F.data(:) = inv( L ) * y;

    case {'pinv'}
      L = [ M ; L ];
      y( end + size(L,1) , end ) = 0;

      F.data(:) = pinv( L )* y ;

    case {'\' , 'backslash'}
      L = [ M ; L ];
      y( end + size(L,1) , end ) = 0;

      F.data(:) = L \ y;
      
    case {'lsqr'}
      L = ( M.'*M + L.'*L );
      y = ( M.' * y );
      
      start_time = clock;

      [F.data(:), flag,relres] = lsqr( L , y , [], 20, [], [], init);
      fprintf('relres: %g\n',relres);
      while flag ~= 0 && flag ~= 3 && flag ~= 4
        [F.data(:), flag,relres] = lsqr( L , y , [], 100, [], [], F.data(:) );
        fprintf('relres: %g\n',relres);
        if etime(start_time) > max_time, flag = -1; break; end
      end
      if flag ~= 0
        warning('it does not converge!!!  (flag: %d)',flag);
      end

    case {'lsqrrec'}
      L = [ M ; L ];
      y( end + size(L,1) , end ) = 0;
      
      start_time = clock;

      [F.data(:), flag,relres] = lsqr( L , y , [], 20, [], [], init);
      fprintf('relres: %g\n',relres);
      while flag ~= 0 && flag ~= 3 && flag ~= 4
        [F.data(:), flag,relres] = lsqr( L , y , [], 100, [], [], F.data(:) );
        fprintf('relres: %g\n',relres);
        if etime(start_time) > max_time, flag = -1; break; end
      end
      if flag ~= 0
        warning('it does not converge!!!  (flag: %d)',flag);
      end
      
  end
  
  
  if NORMALIZE
    F.data(:) = bsxfun( @times , reshape( F.data , prod(size(F,1:3)) , T ) , y_scale );
  end
  F.data(:) = bsxfun( @plus , reshape( F.data , prod(size(F,1:3)) , T ) , y_mean );

  
  if nargout > 1
    
    E(1) = sum( ( Mo*F.data(:) - yo ).^2 ) / N;
    E(2) = sum( ( Lo*F.data(:)      ).^2 ) * dV;
    
  end


end
