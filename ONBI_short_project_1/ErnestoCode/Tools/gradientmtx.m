function G = gradientmtx( sz , dim , varargin )


if 0
%%

stencil =  'quadratic';
bounds  = {'decay',5};

x = 20+[ 0 5 5.5 6 7 8 10 10.5 11 11.5 12 12.5 15 16 16.5 18.5 20];
y = [ cumsum( [ones(1,numel(x)-6)*2  zeros(1,2) -ones(1,3) ] ) 10 ];
y(1)=0;
xx = linspace( x(1)-10 , x(end)+10, 1000 );

cla; hold on
plot( x , y , '*r' ); 
plot( xx,Interp1D(y(:),x,xx,bounds{:},'outside_value',nan) ,'-g' );
 

set(gca,'xtick',dualVector(x),'xgrid','on'); 

G = gradientmtx( numel(x) , 1 , x , stencil , bounds{:} );
d = G*y(:);

plot(  vec( [ x - 0.3; x + 0.3 ; x*0 + NaN ] ) , vec( [ y - d'*0.3 ; y + d'*0.3  ; y ] ) , '-m' ,'linewidth',2 )

hold off; axis equal 

%%

end


%{

[xx,yy,zz,tt] = ndgrid( 1:1:20 , 1:1:25 , 1:1:35 , 1:5 );

v  = xx + yy.^2 + 4*zz + tt.^4;
vx = xx*0 + 1;
vy = 2*yy;
vz = 4*zz*0 + 4;
vt = 4*tt.^3;

%}

  if nargin < 2
    dim = 1;
  end
  if ischar( dim ) || iscell( dim )
    varargin{end+1} = dim;
    dim = 1;
  end

  
  METHOD = 'centered';
  BOUNDS = 'none';
  
  
  [varargin,METHOD] = parseargs(varargin,'Forward'          ,'$FORCE$',{'forward'         ,METHOD});
  [varargin,METHOD] = parseargs(varargin,'Backward'         ,'$FORCE$',{'backward'        ,METHOD});
  [varargin,METHOD] = parseargs(varargin,'Centered'         ,'$FORCE$',{'centered'        ,METHOD});
  [varargin,METHOD] = parseargs(varargin,'SUBDifferential','Subdiff'  ,'$FORCE$',{'subdifferential' ,METHOD});
  [varargin,METHOD] = parseargs(varargin,'Quadratic'        ,'$FORCE$',{'quadratic'       ,METHOD});
  
  
  [varargin,BOUNDS] = parseargs(varargin,'none','value' ,'$FORCE$',{'none'       ,BOUNDS});
  [varargin,BOUNDS] = parseargs(varargin,'periodic'     ,'$FORCE$',{'periodic'   ,BOUNDS});
  [varargin,BOUNDS] = parseargs(varargin,'symmetric'    ,'$FORCE$',{'symmetric'  ,BOUNDS});
  [varargin,BOUNDS] = parseargs(varargin,'closest'      ,'$FORCE$',{'closest'    ,BOUNDS});

  [varargin,DEC,dec_distance] = parseargs(varargin,'decay','$DEFS$',1);
  if DEC
    BOUNDS = 'decay';
  end
  
  
  d  = sz( dim );
  
  if numel( varargin )
    
    if isnumeric( varargin{1} )
      X = varargin{1};
    elseif iscell( varargin{1} )
      X = ( 1:d ) * varargin{1}{1};
    end
    varargin(1) = [];

  else

    X = 1:d;

  end
  X = X(:).';
  
  if numel(X) ~= d
    error('incorrect coordinates');
  end
  

  if d == 1
    
    G1 = sparse( 0 );
    G = kron( kron( eye(prod(sz(dim+1:end))) , G1 ) , eye(prod(sz(1:dim-1))) );
    
  elseif d > 1

    switch lower(BOUNDS)
      case 'periodic'
        W = X - ( X(d) - X(d-1) )/2 - ( X(2) - X(1) )/2 - X(d) + X(1);
        Y = X + ( X(d) - X(d-1) )/2 + ( X(2) - X(1) )/2 + X(d) - X(1);
      case 'symmetric'
        W = [0 cumsum(fliplr(diff(X)))];
        W = X(1) - (X(2)-X(1)) - W(d) + W;
        Y = X(d) + X(d)-X(d-1) + [0 cumsum(fliplr(diff(X)))];
      case 'decay'
        W = linspace( -( X(d) - X(1) ) , 0 , d ) + X( 1 ) - dec_distance;
        Y = linspace( 0 , ( X(d) - X(1) )  , d ) + X(end) + dec_distance;
    end
    
    switch lower( METHOD )
      case {'forward'}
        G1 = sparse( 1:d-1 , 1:d-1 , -1 ./( X(2:d) - X(1:d-1) ) , d , d , 2*d ) + ...
             sparse( 1:d-1 , 2:d   ,  1 ./( X(2:d) - X(1:d-1) ) , d , d , 2*d );
        switch lower(BOUNDS)
          case {'none'}
            G1( d , [ d  d-1] )= [1 -1]./( X( d ) - X(d-1) );
          case {'periodic'}
            G1( d , [ 1   d ] )= [1 -1]./( Y( 1 ) - X( d ) );
          case {'symmetric'}
            %G1( d , [ d   d ] )= [1 -1]./( Y( 1 ) - X( d ) );
          case {'closest'}
            %G1( d , [ d   d ] )= [1 -1]./( Y( 1 ) - X( d ) );
          case {'decay'}
            G1( d , [ d     ] )= [  -1]./( Y( 1 ) - X( d ) );
        end

      case {'backward'}
        G1 = sparse( 2:d , 1:d-1 , -1 ./( X(2:d) - X(1:d-1) ) , d , d , 2*d ) + ...
             sparse( 2:d , 2:d   ,  1 ./( X(2:d) - X(1:d-1) ) , d , d , 2*d );
        switch lower(BOUNDS)
          case {'none'}
            G1( 1 , [ 2   1 ] )= [1 -1]./( X( 2 ) - X( 1 ) );
          case {'periodic'}
            G1( 1 , [ 1   d ] )= [1 -1]./( X( 1 ) - W( d ) );
          case {'symmetric'}
          case {'closest'}
          case {'decay' }
            G1( 1 , [ 1     ] )= [1   ]./( X( 1 ) - W( d ) );
            
        end

      case {'centered'}
        G1 = sparse( 2:d-1 , 1:d-2 , -1 ./( X(3:d) - X(1:d-2) ) , d , d , 2*d ) + ...
             sparse( 2:d-1 , 3:d   ,  1 ./( X(3:d) - X(1:d-2) ) , d , d , 2*d );
        switch lower(BOUNDS)
          case {'none'}
            G1(  1  , [ 2   1 ] )= [1 -1]./( X( 2 ) - X( 1 ) );
            G1(  d  , [ d  d-1] )= [1 -1]./( X( d ) - X(d-1) );
          case {'periodic'}
            G1(  1  , [ 2   d ] )= [1 -1]./( X( 2 ) - W( d ) );
            G1(  d  , [ 1  d-1] )= [1 -1]./( Y( 1 ) - X(d-1) );
          case {'symmetric'}
            G1(  1  , [ 2   1 ] )= [1 -1]./( X( 2 ) - W( d ) );
            G1(  d  , [ d  d-1] )= [1 -1]./( Y( 1 ) - X(d-1) );
          case {'closest'}
            G1(  1  , [ 2   1 ] )= [1 -1]./( X( 2 ) - X( 1 ) )/2;
            G1(  d  , [ d  d-1] )= [1 -1]./( X( d ) - X(d-1) )/2;
          case {'decay'}
            G1(  1  , [ 2     ] )= [1   ]./( X( 2 ) - W( d ) );
            G1(  d  , [    d-1] )= [  -1]./( Y( 1 ) - X(d-1) );
        end


      case {'subdifferential'}
        G1 = ( gradientmtx( d , 1 , X , 'forward'  , BOUNDS , dec_distance ) + ...
               gradientmtx( d , 1 , X , 'backward' , BOUNDS , dec_distance ) )/2;

             
             
      case {'quadratic'}
        G1 = sparse( 2:d-1 , 1:d-2 , ( 1./( X(3: d )-X(1:d-2) ) + 1./( X(1:d-2) - X(2:d-1) ) ) , d , d , 2*d ) + ...
             sparse( 2:d-1 , 2:d-1 , ( 1./( X(2:d-1)-X(1:d-2) ) + 1./( X(2:d-1) - X(3: d ) ) ) , d , d , 2*d ) + ...
             sparse( 2:d-1 , 3: d  , ( 1./( X(1:d-2)-X(3: d ) ) + 1./( X(3: d ) - X(2:d-1) ) ) , d , d , 2*d ) ;
      
        switch lower(BOUNDS)
          case {'none'}
            G1(  1  , 1 ) =  1/(X( 1 )-X( 2 )) + 1/(X( 1 )-X( 3 ));
            G1(  1  , 2 ) = -1/(X( 1 )-X( 2 )) - 1/(X( 2 )-X( 3 ));
            G1(  1  , 3 ) = -1/(X( 1 )-X( 3 )) + 1/(X( 2 )-X( 3 ));

            G1(  d  ,d-2) = -1/(X(d-2)-X(d-1)) + 1/(X(d-2)-X( d ));
            G1(  d  ,d-1) =  1/(X(d-2)-X(d-1)) + 1/(X(d-1)-X( d ));
            G1(  d  , d ) = -1/(X(d-2)-X( d )) - 1/(X(d-1)-X( d ));
          case {'periodic'}
            G1(  1  , 1 ) = -1/(W( d )-X( 1 )) + 1/(X( 1 )-X( 2 ));
            G1(  1  , 2 ) =  1/(W( d )-X( 2 )) - 1/(X( 1 )-X( 2 ));
            G1(  1  , d ) =  1/(W( d )-X( 1 )) - 1/(W( d )-X( 2 ));

            G1(  d  ,d-1) = -1/(X( d )-X(d-1)) - 1/(X(d-1)-Y( 1 ));
            G1(  d  , d ) =  1/(X( d )-X(d-1)) + 1/(X( d )-Y( 1 ));
            G1(  d  , 1 ) = -1/(X( d )-Y( 1 )) + 1/(X(d-1)-Y( 1 ));
          case {'symmetric'}
            G1(  1  , 1 ) =  1/(X( 1 )-X( 2 )) + 1/(X( 2 )-W( d ));
            G1(  1  , 2 ) =  1/(X( 2 )-X( 1 )) + 1/(W( d )-X( 2 ));

            G1(  d  ,d-1) =  1/(X(d-1)-X( d )) + 1/(Y( 1 )-X(d-1));
            G1(  d  , d ) =  1/(X( d )-X(d-1)) + 1/(X(d-1)-Y( 1 ));
          case {'closest'}
            G1(  1  , 1 ) =  1/(X( 1 )-X( 2 )) + 1/(X( 2 )-W( d ));
            G1(  1  , 2 ) =  1/(X( 2 )-X( 1 )) + 1/(W( d )-X( 2 ));

            G1(  d  ,d-1) =  1/(X(d-1)-X( d )) + 1/(Y( 1 )-X(d-1));
            G1(  d  , d ) =  1/(X( d )-X(d-1)) + 1/(X(d-1)-Y( 1 ));
          case {'decay'}
            G1(  1  , 1 ) =  1/(X( 1 )-W( d )) + 1/(X( 1 )-X( 2 ));
            G1(  1  , 2 ) =  1/(W( d )-X( 2 )) - 1/(X( 1 )-X( 2 ));

            G1(  d  ,d-1) = -1/(X( d )-X(d-1)) + 1/(Y( 1 )-X(d-1));
            G1(  d  , d ) =  1/(X( d )-X(d-1)) + 1/(X( d )-Y( 1 ));
        end

    end
    G = kron( kron( eye(prod(sz(dim+1:end))) , G1 ) , eye(prod(sz(1:dim-1))) );

  else    
    
    G = @(v,x) cat(1 ,                                                     ...
                   (                                                       ...
                     v(2:end  ,:,:,:,:,:,:,:,:,:,:,:,:,:,:) -              ...
                     v(1:end-1,:,:,:,:,:,:,:,:,:,:,:,:,:,:)                ...
                   ) ./  ( x(2:end) - x(1:end-1) )                         ...
                   ,                                                       ...
                   (                                                       ...
                     v(  end  ,:,:,:,:,:,:,:,:,:,:,:,:,:,:) -              ...
                     v(  end-1,:,:,:,:,:,:,:,:,:,:,:,:,:,:)                ...
                   ) ./  ( x(end) - x(end-1) )                             ...
                  );
    
  end











%   nd = numel( sz );
%   sz = [ sz(:)' 1 1 ];
%   sz( max( [ numel(varargin)  nargout+1   numel(sz) ] ) ) = 1;
%   sz( ~sz ) = 0;
%   N  = prod( sz );
% 
%   
%   if numel( varargin ) == 1
% 
%     H = varargin{1};
%     if isscalar( H )
%       H = repmat(H,1,nd);
%     end
%     if numel(H) < nd
%       H(nd) = 1; H( ~H ) = 1;
%     end
%     Xs = arrayfun( @(n) ( 1:sz(n) )*H( n ) , 1:nd , 'un' , false );
%     
%   else
% 
%     Xs = varargin;
%     for d = find( cellfun(@isempty,Xs) )
%       Xs{d} = 1:sz(d);
%     end
%     for d = numel(Xs)+1:nargout
%       Xs{d} = 1:sz(d);
%     end
%     
%     sizeXs = cellfun( @(x) numel(x) , Xs );
%     if any( sz(1:numel(sizeXs)) - sizeXs )
%       error('gradientmtx:InconsistentSizes','Numbers of coordinates have to equal sizes');
%     end
%     
%   end
%   Xs = cellfun( @(x) x(:) , Xs , 'un' , false );
%   
%   
%   idxs = ndmat( num2cell( sz ) , 'cell' , 'nocat' );
% %   bytessize('idxs')
% %   idxs = cellfun( @(id) uint16(id) , idxs , 'UniformOutput', false );
% %   bytessize('idxs')
%   idxs = cellfun( @(id) id(:) , idxs , 'UniformOutput', false );
%   
%   for d = 1:nargout
%     if d > nd
%       varargout{d} = sparse( N , N );
%       continue;
%     end
%     
%     M = sparse( N , N );
%     
%     r = sz;  r(d) = 1;
%     
%     %forward
%     idxs_ = idxs;
%     idxs_{d} = idxs_{d} + 1;
%     idxs_{d}( idxs_{d} > sz(d) ) = sz(d);
% 
%     deltas = 1./diff( Xs{d} )/2;
%     deltas = ipermute( [ deltas(1)*2 ; deltas(2:end) ; deltas(end)*2 ] , [ d 1:d-1 d+1:nd ] );
%     deltas = repmat( deltas , r );
%     
%     M = M + sparse( 1:N , sub2ind(sz,idxs_{:}) , deltas(:) ,N,N );
%     M = M - sparse( 1:N ,     1:N              , deltas(:) ,N,N );
% 
%     
%     %backward
%     idxs_ = idxs;
%     idxs_{d} = idxs_{d} - 1;
%     idxs_{d}( idxs_{d} < 1 ) = 1;
% 
%     deltas = 1./diff( Xs{d} )/2;
%     deltas = ipermute( [ deltas(1)*2 ; deltas(1:end-1) ; deltas(end)*2 ] , [ d 1:d-1 d+1:nd ] );
%     deltas = repmat( deltas , r );
% 
%     M = M - sparse( 1:N , sub2ind(sz,idxs_{:}) , deltas(:) ,N,N );
%     M = M + sparse( 1:N ,    1:N               , deltas(:) ,N,N );
% 
%     varargout{d} = M;
%   end

end
