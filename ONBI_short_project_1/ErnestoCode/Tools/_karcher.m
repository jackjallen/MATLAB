function M = karcher( X , EXP , LOG , varargin )

  
  M = X{1};

  while 1
    Mp = M;
    
    LOGS = cellfun( @(xi) vec( LOG(M,xi) ) , X , 'UniformOutput',false);
    M = EXP( M , mean( cell2mat( LOGS ) , 2 ) );
    
%     maxnorm( Mp - M )
    if maxnorm( Mp - M ) < 1e-15, break; end
  end
  
  






























%   if ischar( EXP )
%     switch EXP
%       case 'DTI'
%         EXP= @(I,x) sqrtm(I)*expm( inv(sqrtm(I))*x*inv(sqrtm(I)) )*sqrtm(I);
%       case 'S'
%         EXP= @(I,x) expS(I,x);
%     end
%   end
%   
%   if ischar( LOG )
%     switch LOG
%       case 'DTI'
%         LOG= @(I,x) sqrtm(I)*logm( inv(sqrtm(I))*x*inv(sqrtm(I)) )*sqrtm(I);
%       case 'S'
%         LOG= @(I,x) logS(I,x);
%     end
%   end
% 
%   N= numel(x);
%   [varargin,i,weigths]= parseargs( varargin , 'Weigths','$DEFS$',ones(N,1) );
%   
%   M     = x{1};
%   yprev = x{1}+1000;
%   
%   iter= 0;
%   while norm( M-yprev ) > 1e-10
%     iter= iter+1;
%     yprev = M;
%     
%     sumlog = 0;
%     for i=1:N
%       sumlog = sumlog + weigths(i)*feval( LOG , M , x{i} );
%     end
%     sumlog = sumlog/ sum( weigths(1:N) );
%     
%     M = feval( EXP , M , sumlog );
% %     norm( M-yprev )
%   end
% %   disp(iter);
% 
% 
%   function x= logS(I,x)
%     Nb = null(I');
% %     t  = acos( I'*x );
%     if I'*x < 0
%       t = pi - 2*asin( norm(I+x)/2 );
%     else
%       t = 2*asin( norm(I-x)/2 );
%     end
%     if t
%       x  = Nb'*x;
%       x  = x/norm(x)*t;
%       x  = Nb*x + I;
%     else
%       x  = I;
%     end
%   end
%   function x= expS(I,x)
%     t = norm( I-x );
%     x = x-I;
%     x = x/norm(x);
%     x = I*cos(t) + sin(t)*x;
%   end

end
