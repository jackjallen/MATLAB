function [L,e] = reallogm( M , L )

  if nargin < 2
    setWarning('off','MATLAB:funm:nonPosRealEig');
    L = logm( M );
    restoreWarning(  'MATLAB:funm:nonPosRealEig');

    if isreal( L )
      if nargout > 1, e = max( abs( vec( M - expm(L) ) ) ); end
      return;
    end

    if max( abs( vec( imag( L ) ) ) ) < 1e-16
      L = real( L );
      if nargout > 1, e = max( abs( vec( M - expm(L) ) ) ); end
      return;
    end




    T = zeros( size(M) );

    indxs_p = ~~triu( ones(size(M)) ,  1 );
    indxs_n = ~~tril( ones(size(M)) , -1 );
    try
      d_mT_t = NumericalDiff( @(t) mT(t) , zeros(sum(indxs_p(:)),1) , 'i' );
    end





    if size(M,1) == 3

      R = M/abs( power(det(M),1/3) );

      try
        q = mat2quat( R );
        nq = min( norm( q(2:4) ) , 1 );
        L = q(2:4)/nq*asin(nq)*2;
        L  = [0 , -L(3) , +L(2); L(3) , 0 , -L(1); -L(2) , L(1) , 0];

      catch

        ry = asin( R(1,3) );
        if ry < pi/2
          if ry > -pi/2
            rx = atan2( -R(2,3) , R(3,3) );
            rz = atan2( -R(1,2) , R(1,1) );
          else
            rx = -atan2( R(2,1) , R(2,2) );
            rz = 0;
          end
        else
          rx = atan2( R(2,1) , R(2,2) );
          rz = 0;
        end

        rx = real(rx);
        ry = real(ry);
        rz = real(rz);

        s = 1e-8;
        while 1
          R = maketransform('radians','rz',rz-s,'ry',ry-s,'rx',rx-s);
          R = R(1:3,1:3);

          setWarning('off','MATLAB:funm:nonPosRealEig');
          L = logm( R );
          restoreWarning(  'MATLAB:funm:nonPosRealEig');

          if isreal( L )
            break;
          end
          s = s*2;
        end

      end

      L_init = L([4 7 8]);

    else

      L_init = zeros(sum(indxs_p(:)),1);

    end

    try
      L = Optimize( @(t) error_rot(t) , L_init , 'methods' , {'conjugate',50,'quasinewton',50,'coordinate',1} , 'ls' , {'quadratic','golden'} ,'verbose',0,'noplot');

      d = L/norm(L)*( 2*pi );

        LL = L; it = 0;
        while norm( LL ) > pi && it < 5
          it = it+1;
          LL = LL - d; 
        end
      if norm( LL ) > pi
        LL = L; it = 0;
        while norm( LL ) > pi && it < 5
          it = it+1;
          LL = LL + d; 
        end
      end
      if norm( LL ) > pi
        LL = L;
      end


      L = LL;

      L = mT(L);
    end
    
  else
    
    if ~isequal( size(M) , size(L) )
      error('bad initialization');
    end
    
  end

  try
    L = Optimize( @(t) error_full(t) , L , 'methods' , {'conjugate',50,'quasinewton',50,'coordinate',1} ,'ls',{'quadratic','golden'} ,'verbose',0,'plot');  
  end
  
  
  
  
  if nargout > 1, e = max( abs( vec( M - expm(L) ) ) ); end
  
  %vfunction 'TT = mT(t), T( indxs_p ) =  t(:); T( indxs_n ) = -t(:); TT = T; end'


  
  
  function [E,J] = error_rot( L )
    mT_L = mT(L);
    e_mT_L = expm( mT_L );
    D = e_mT_L - M;
    E = D(:).'*D(:);
    
    if nargout > 1
      J = 2*D(:).'*d_expm( mT_L )*d_mT_t;
    end
  end
  
  function [E,J] = error_full( L )
    e_L = expm( L );
    D = e_L - M;
    E = D(:).'*D(:);

    if nargout > 1
      J = 2*D(:).'*d_expm( L );
    end
  end



    function TT = mT(t)
      T( indxs_p ) =  t(:);
      T( indxs_n ) = -t(:);
      TT = T;
    end

end
