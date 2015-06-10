function [ X , U , R ] = shoot( x0 , u0 , G , FLOW , varargin )

  if nargin < 4 || isempty( FLOW ), FLOW = 1; end

  dG = [];
  if iscell( G ) && numel( G ) > 1
    dG = G{2};
    G  = G{1};
  elseif iscell(G)
    G  = G{1};
  end
  
  if ~isa( G , 'function_handle' )
    error('G is expected to be a function_handle');
  end
  
  if isempty( dG )
    try
      [g,dg] = G(x0(:));
    catch
      dG = @(x) NumericalDiff( @(z) G(z) , x , 'c' );
    end
  else
    if ~isa( dG , 'function_handle' )
      error('dG is expected to be a function_handle');
    end
  end

  
  
  if numel(x0) ~= numel(u0), error('initial pos and initial vel must have equal sizes'); end
  
  N = numel(x0);
  
  x0 = x0(:);
  u0 = u0(:);
  
  
  
  FLOW = [0,FLOW(:)'];
  if ~issorted( FLOW ), error('times must be ascending'); end

  X = zeros( N , numel(FLOW)-1 );
  if nargout > 1
    U = zeros( N , numel(FLOW)-1 );
  end

  
  TF = max(FLOW);
  
  SHOW_PROGRESS_BAR = true;
  [varargin,SHOW_PROGRESS_BAR] = parseargs( varargin , 'NoProgressBar' , '$FORCE$' , {false , SHOW_PROGRESS_BAR} );
  [varargin,SHOW_PROGRESS_BAR] = parseargs( varargin , 'ProgressBar'   , '$FORCE$' , {true  , SHOW_PROGRESS_BAR} );
  if SHOW_PROGRESS_BAR
    progressBar;
    progressBar( 0 ,'title' );
    CLEANUP = onCleanup( @() progressBar(1) );
  end
  
  
  
  METHOD = 'optimal_control';
  [varargin,METHOD] = parseargs( varargin , 'OptimalControl' , '$FORCE$' , {'optimal_control' , METHOD} );
  [varargin,METHOD] = parseargs( varargin , 'EulerLagrange'  , '$FORCE$' , {'euler_lagrange'  , METHOD} );
  

  ode_opts = odeset( 'RelTol',100*eps(1),'AbsTol',eps(1),'NormControl','off','InitialStep',[],'MaxStep',max(FLOW)/2 ,'Refine',1,'stats','off', varargin{:} );

  if strcmpi( METHOD , 'optimal_control' )
  
    p0 = -( G(x0) * u0 );
    R = struct('t',[],'x',[],'p',[],'u',[],'H',[],'G',[]);
    for i = 1:(numel(FLOW)-1)
      if FLOW(i+1) > FLOW(i)
        [T,xp] = ode45( @(t,xp) xp_dot(xp,t) , [FLOW(i) FLOW(i+1)] , [x0;p0] , ode_opts );

        x0 = xp(end,     1:N   );  x0 = x0(:);
        p0 = xp(end, (N+1):end );  p0 = p0(:);
      else
        T  = FLOW(i);
        xp = [x0;p0].';
      end
      X(:,i) = x0;
      if nargout > 1
        G_x = G( x0 );
        U(:,i) = - ( G_x \ p0 );
        U(:,i) = - pinv( G_x ) * p0;
      end

      if nargout > 2
        id = numel( R.t );

        R.t = [ R.t , T(:).'              ];
        R.x = [ R.x , xp(:,     1:N   ).' ];
        R.p = [ R.p , xp(:, (N+1):end ).' ];
        R.u( N , numel( R.t ) ) = 0;
        R.H( 1 , numel( R.t ) ) = 0;
        R.G( numel(x0) , numel(x0) , numel(R.t) ) = 0;

        for t = (id+1):numel(R.t)
          G_x = G( R.x(:,t) );

          R.G(:,:,t) = G_x;
          %R.u(:,t)   = -( G_x \ R.p(:,t) );
          R.u(:,t)   = - pinv( G_x ) * R.p(:,t);

          R.H(1,t) = ( R.u(:,t).' * G_x * R.u(:,t) )/2  +  R.p(:,t).' * R.u(:,t);
          R.arc_length(1,t) = sqrt(  R.u(:,t).' * G_x * R.u(:,t) );
        end
      end

    end

  elseif strcmpi( METHOD , 'euler_lagrange' )

    
    R = struct('t',[],'x',[],'u',[],'u_dot',[],'G',[]);
    for i = 1:(numel(FLOW)-1)
      if FLOW(i+1) > FLOW(i)
        [T,xu] = ode45( @(t,xu) xu_dot(xu,t) , [FLOW(i) FLOW(i+1)] , [x0;u0] , ode_opts );

        x0 = xu(end,     1:N   );  x0 = x0(:);
        u0 = xu(end, (N+1):end );  u0 = u0(:);
      else
        T  = FLOW(i);
        xu = [x0;u0].';
      end
      X(:,i) = x0;
      if nargout > 1
        U(:,i) = u0;
      end

      if nargout > 2
        id = numel( R.t );

        R.t = [ R.t , T(:).'              ];
        R.x = [ R.x , xu(:,     1:N   ).' ];
        R.u = [ R.u , xu(:, (N+1):end ).' ];
        R.u_dot( N , numel( R.t ) ) = 0;
        R.G( numel(x0) , numel(x0) , numel(R.t) ) = 0;

        for t = (id+1):numel(R.t)
          G_x = G( R.x(:,t) );
          R.G(:,:,t) = G_x;
          
          v = xu_dot( [ R.x(:,t) ; R.u(:,t) ] );
          R.u_dot(:,t) = v( (N+1):end );
        end
      end

    end
    
  else
    error('unknow METHOD');
  end


  
  
  function v = xu_dot( xu , t )
    if nargin > 1  &&  SHOW_PROGRESS_BAR
      progressBar( t/TF , '%T' )
    end
    
    x = xu(     1:N   );  x = x(:);
    u = xu( (N+1):end );  u = u(:);
    
    if isempty( dG )
      [ G_x , dG_x ] = G( x );
    else
       G_x =  G( x );
      dG_x = dG( x );
    end
    
    u_dot = ( dG_x.' * kron( u , u )/2 - kron( u.' , eye(N) ) * dG_x * u );
    u_dot = G_x \ u_dot;
    
    v = [ u ; u_dot ];
  end


  function v = xp_dot( xp , t )
    if nargin > 1  &&  SHOW_PROGRESS_BAR
      progressBar( t/TF , '%T' )
    end
    
    x = xp(     1:N   );  x = x(:);
    p = xp( (N+1):end );  p = p(:);
    
    if isempty( dG )
      [ G_x , dG_x ] = G( x );
    else
       G_x =  G( x );
      dG_x = dG( x );
    end
    
    %u = -( G_x \ p );  u = u(:);
    u = - pinv( G_x ) * p;  u = u(:);
    
    v = [ u ; - dG_x.' * vec( u * u.' )/2 ];
  end

end
