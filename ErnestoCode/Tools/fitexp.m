function [ lambda , alpha , beta ] = fitexp( T , Y , add_constant , varargin )
% 
%   Y ~ alpha * exp( -lambda * T ) + beta
% 
%
%
% [ lambda , alpha , beta ] = fitexp( T , Y , add_constant(true or false) , weights , cost_fun , tune )
%
%

  if nargin < 3, add_constant = true;  end

  if ~isvector( T ), error( 'a vector expected as first  arg' ); end
  if ~isvector( Y ), error( 'a vector expected as second arg' ); end
  
  T = T(:);
  Y = Y(:);
  
  
  if numel( T ) ~= numel( Y )
    error('different sizes');
  end
  
  LIMITS = maxnum('single');
  R = range( T );
  if R(1) < 0
    LIMITS( end+1 ) = -709.7/R(1);   %%avoid underflow: exp(  709 ) = Inf
  end
  if R(1) > 0
    LIMITS( end+1 ) =  746/R(1);   %%avoid underflow: exp( -746 ) = 0
  end

  LIMITS = [ 0 , min( LIMITS ) ];
  
  
  if add_constant
    lambda = Optimize1D( @(g) ewlr( { exp(-g*T) , 1 } , Y , varargin{:} ) , LIMITS , [] , 'METHODS' , 'exhaustive' , 'exhaustive.N',1000);

    if nargout > 1
      ab = wlr( { exp(-lambda*T) , 1 } , Y , varargin{:} );
      alpha = ab(1); beta  = ab(2);
    end
    
  else
    lambda = Optimize1D( @(g) ewlr(   exp(-g*T)       , Y , varargin{:} ) , LIMITS , [] , 'METHODS' , 'exhaustive' , 'exhaustive.N',1000);

    if nargout > 1
      alpha = wlr( exp(-lambda*T) , Y , varargin{:} );
      beta  = 0;
    end

  end
  
  
%   GS = [0,geospace(1e-20,LIMITS(2),2000)];
%   figure;
%   subplot(2,1,1);
%   for i = 1:numel(GS)
%     [JJ(i),GG(i)] = expfit_0( GS(i) );
%   end
%   subplot(2,1,1);  semilogx( GS , JJ , '.-r' );
%   subplot(2,1,2);  semilogx( GS , GG , '.-r' );
%   subplot(2,1,1);  vline( lambda );
%   subplot(2,1,2);  vline( lambda ); hline( 0 );
  
  
  
  
  
%   function [J,dJ] = expfit_0( L )
%     
%     Z = exp( -L * T );
%     alpha_ = ( Y.' * Z ) / ( Z.' * Z );
%     
%     ZZ = alpha_ * Z - Y;
%     
%     J = sum( ZZ .^ 2 );
%     
%     %JJ = @(L) sum(  ( ( Y.' * exp( -L * T ) ) / ( exp( -L * T ).' * exp( -L * T ) ) * exp( -L * T ) - Y ).^2 );
%     %GG = @(L) -2*( Y.' * exp(-L*T) )/( exp(-L*T).' * exp(-L*T) )*( ( Y.' * exp(-L*T) )/( exp(-L*T).' * exp(-L*T) ) * exp(-L*T) - Y ).'*( exp(-L*T) .* T );
% 
%     if nargout > 0
%       dJ = - 2 * alpha_ * ZZ.' * ( T .* Z );
%     end
%     
%   end


%   function [J,dJ] = expfit_1( L )
%     
%     Z = exp( -L * T );
%     
%     [ab,J] = wlr( {Z,1} , Y );
%     
%     if nargout > 1
%       dJ = -2 * ab(1) * ( ab(1) * Z + ab(2) - Y ).' * ( T .* Z );
%     end
% 
%   end
    



end
