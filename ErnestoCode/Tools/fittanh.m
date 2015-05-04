function [F,ener] = fittanh( x , y , bounds )
if 0
  
  x = randn(10000,1); x = x(:);
  
  A = 0.5; B = 20; C = 1; D = 0;
  x = x + B*0.75;
  
  y = C * tanh( A*( x-B ) ) + D;
  y = y + randn(size(y))/20 + 0.4;
  
  F = fittanh( x , y , [-1 1] );
  
  set( ezplot( F , centerscale( range(x) , 4 ) ) , 'color','g');
  hold on;
  set( ezplot( @(x) C*tanh( A*(x-B) ) + D , centerscale( range(x) , 4 ) ) , 'color','m')
  hplot( x , y , '.r' )
  hold off
  
end


  if nargin < 3
    bounds = [-inf,inf];
    bounds = [-1,1];
  end

  if numel( x ) ~= numel( y ), error('x and y must have the same number of points'); end
  if numel(bounds) ~= 2,       error('bounds is expected to be 2 elements'); end
  if bounds(1) >= bounds(2),   error('bounds must be sorted'); end
  
  w = ~isnan(x) & ~isnan(y) & ~isinf(x) & ~isinf(y);
  
  x = x(w);
  y = y(w);
  
  if numel(x) < 3, error('too few points'); end
  

  
  ENERGY = @(abcd) fro2( abcd(3)*tanh( exp(abcd(1))*( x - abcd(2) ) ) + abcd(4) - y );
  
  if ~isinf( bounds(1) ) && ~isinf( bounds(2) )

    C = ( bounds(2)-bounds(1) )/2;
    D = ( bounds(2)+bounds(1) )/2;
    ABCD = [ 0 0 C D ];
    
    
    ABCD(1:2) = ExhaustiveSearch( @(ab) ENERGY([ab,C,D]) , [0 0] , 20 , 45 , 'maxTIME',10*60,'verbose');
    

    
%     
%     ener = NaN;
%     
%     while true
%       ener_old = ener;
%       
%       ABCD(2) = Optimize1D( @(b) fro2( ABCD(3)*tanh( exp(ABCD(1)) * ( x - b       ) ) + ABCD(4) - y ) , [ -2 2 ]*max(ABCD(2),1) , [] , 'METHODS' , 'exhaustive' , 'exhaustive.N',11 , 'exhaustive.MAX_ITS',61);
%       ABCD(1) = Optimize1D( @(a) fro2( ABCD(3)*tanh( exp(a)       * ( x - ABCD(2) ) ) + ABCD(4) - y ) , [ -2 2 ]*max(ABCD(1),1) , [] , 'METHODS' , 'exhaustive' , 'exhaustive.N',11 , 'exhaustive.MAX_ITS',61);
%       
%       ener = fro2( ABCD(3)*tanh( exp(ABCD(1))   * ( x - ABCD(2) ) ) + ABCD(4) - y );
%       
% %       %{
%       figure(1000);
%       ezplot( @(x)ABCD(3)*tanh(exp(ABCD(1))*(x-ABCD(2)))+ABCD(4) , [ centerscale( range(x) , 4 ) , range(y) ] ) ; hplot(x,y,'.r')
%       title( uneval(  ener ) );
%       drawnow;
% %       %}
%       
%       if ener_old <= ener
%         break;
%       end
%     end
    

  else
    
    
  end


  F = eval( sprintf( '@(x) %s * tanh( %s * ( x - %s ) ) + %s' , uneval(ABCD(3)) , uneval(exp(ABCD(1))) , uneval(ABCD(2)) , uneval(ABCD(4)) ) );
  

  
  function [e,d] = ener1( A , B , C , D )
    
    Ax_B = A * ( x - B );
    
    fx = C * tanh( Ax_B ) + D;
    
    r  = fx - y;
    
    e  = r.' * r;
    
    if nargout > 1
      
      d = [0 0];
      
      sech_Ax_B_2_r = sech( Ax_B ).^2 .* r;
      d(1) = 2 * C * sum( ( x - B ) .* sech_Ax_B_2_r , 1 );
      d(2) = -2 * A * C * sum( sech_Ax_B_2_r , 1 );
      
    end
    
  end
  
  
end
