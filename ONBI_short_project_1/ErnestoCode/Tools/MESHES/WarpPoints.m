function points= WarpPoints( points , o , t , varargin)
% 
% points= WarpPoints( points, origin_points , target_points ,
%                               'NumberSpatialDimensions', [2 3]
%                               'Plot' 
%                               'BASIS', [ 'R','RLOGR' ]
% 

  nsd= size(o,2);
  plotear= 'no';
  basis = 'r';

  zz= parseargs( varargin, 'NumberSpatialDimensions' ); if zz, nsd= varargin{zz+1}; end
  zz= parseargs( varargin,'Basis' );                         
      if zz
        if      strcmp( lower( varargin{zz+1} ), 'r' )
          basis= 'r';
        elseif  strcmp( lower( varargin{zz+1} ), 'rlogr' )
          basis= 'rlogr';
        end
      end
  zz= parseargs( varargin,'Plot'  );                      if zz, plotear='yes';      end

  R2= R2matrix(o,o);
  if      strcmp( basis, 'rlogr')
    warning('off', 'MATLAB:log:logOfZero');
    K= R2.*log(R2);
    warning('on' , 'MATLAB:log:logOfZero');
    K(isnan(K))=0;
  elseif  strcmp( basis, 'r')
    K= sqrt(R2);
  end

  NPo= size(o,1);
  L= [ K o ones(NPo,1) ; [ o  ones(NPo,1) ]' zeros(nsd+1) ];
  W= pinv(L)*[ t ; zeros(nsd+1,nsd) ];

  points= mexwarp_new( points , o , W , basis);

  if strcmp( plotear , 'yes' )
%   if 1
    figure;
    hold on

    if nsd == 2
        plot( o(:,1),o(:,2), '.g', 'MarkerSize',25 , 'DisplayName','Origin')
        plot( t(:,1),t(:,2), '.m', 'MarkerSize',25 , 'DisplayName','Target')
        onew = mexwarp_new( o , o , W , basis);

        plot( onew(:,1),onew(:,2), '+k', 'MarkerSize',10 , 'DisplayName','Final')
        legend toggle

        for p=1:size( o,1 )
          plot( [ o(p,1) ; onew(p,1) ] , [ o(p,2) ; onew(p,2) ] , '-g' )
        end

        [x,y]= meshgrid( linspace( min(o(:,1)), max(o(:,1)),51 ) , linspace( min(o(:,2)), max(o(:,2)),51 ) );
        p    = [ x(:) y(:) ];
        pnew = mexwarp_new( p , o , W , basis );

        x = reshape( p(:,1), size(x) );
        y = reshape( p(:,2), size(y) );
        xp= reshape( pnew(:,1), size(x) );
        yp= reshape( pnew(:,2), size(y) );
            plot(  x(:,1:5:end)  ,  y(:,1:5:end)  , ':g' )
            plot(  x(1:5:end,:)' ,  y(1:5:end,:)' , ':g' )
            plot( xp(:,1:5:end)  , yp(:,1:5:end)  , ':m' )
            plot( xp(1:5:end,:)' , yp(1:5:end,:)' , ':m' )
    else
        plot3( o(:,1),o(:,2),o(:,3), '.g', 'MarkerSize',25 , 'DisplayName','Origin')
        plot3( t(:,1),t(:,2),t(:,3), '.m', 'MarkerSize',25 , 'DisplayName','Target')
        onew = mexwarp_new( o , o , W , basis );

        plot3( onew(:,1),onew(:,2),onew(:,3), '+k', 'MarkerSize',10 , 'DisplayName','Final')
        legend toggle    

        for p=1:size( o,1 )
          plot3( [ o(p,1) ; onew(p,1) ] , [ o(p,2) ; onew(p,2) ] , [ o(p,3) ; onew(p,3) ] , '-g' )
        end

        [x,y,z]= ndgrid( linspace( min(o(:,1)), max(o(:,1)),21 ) , linspace( min(o(:,2)), max(o(:,2)),21 ) , linspace( min(o(:,3)), max(o(:,3)),21 ));
        p    = [ x(:) y(:) z(:) ];

        pnew = mexwarp_new( p , o , W , basis);

        x = reshape( p(:,1), size(x) );
        y = reshape( p(:,2), size(y) );
        z = reshape( p(:,3), size(z) );
        xp= reshape( pnew(:,1), size(x) );
        yp= reshape( pnew(:,2), size(y) );
        zp= reshape( pnew(:,3), size(y) );
        for i=1:5:size( z,3 )
            plot3(  x(:,1:5:end,i)  ,  y(:,1:5:end,i)  ,  z(:,1:5:end,i)  , ':g')
            plot3(  x(1:5:end,:,i)' ,  y(1:5:end,:,i)' ,  z(1:5:end,:,i)' , ':g')
            plot3( xp(:,1:5:end,i)  , yp(:,1:5:end,i)  , zp(:,1:5:end,i)  , ':m')
            plot3( xp(1:5:end,:,i)' , yp(1:5:end,:,i)' , zp(1:5:end,:,i)' , ':m')
        end

        view(3)
        xlabel('X')
        ylabel('Y')
        zlabel('Z')
    end

    hold off
    axis image
  end

function R2= R2matrix( a,b )
  NPa= size(a,1);
  NPb= size(b,1);
  [D1, D2]= meshgrid( a(:,1) , b(:,1) );
  Rx = reshape( [D1;D2] , NPb, 2*NPa )*kron( eye(NPa), [1 ; -1] ) ;
  [D1, D2]= meshgrid( a(:,2) , b(:,2) );
  Ry = reshape( [D1;D2] , NPb, 2*NPa )*kron( eye(NPa), [1 ; -1] ) ;
  if size(a,2) > 2 & size(b,2) > 2
    [D1, D2]= meshgrid( a(:,3) , b(:,3) );
    Rz = reshape( [D1;D2] , NPb, 2*NPa )*kron( eye(NPa), [1 ; -1] ) ;
  else
    Rz=0;
  end
  R2= Rx.^2+Ry.^2+Rz.^2;
