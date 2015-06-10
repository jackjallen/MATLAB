function [ lambda , alpha , beta ] = fitexp_anchor( T , Y , TA , YA )
% 
%   Y ~ alpha * exp( -lambda * T ) + beta
% 


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
  

  lambda = Optimize1D( @(g) ewlr( exp(-g*T) - exp(-g*TA) , Y - YA , 1 ,'sse' ) , LIMITS , [] , 'METHODS' , 'exhaustive' , 'N',1000);

  if nargout > 1
    alpha = wlr( exp(-lambda*T) - exp(-lambda*TA) , Y - YA , 1 , 'sse' );
    beta  = YA - alpha*exp( -lambda*TA );
  end

  
%   GS = [0,geospace(1e-20,LIMITS(2),2000)];
%   figure;
%   for i = 1:numel(GS)
%     JJ(i) = ewlr( exp(-GS(i)*T) - exp(-(i)*TA) , Y - YA , 1 ,'sse' );
%   end
%   semilogx( GS , JJ , '.-r' );
%   vline( lambda );

end


