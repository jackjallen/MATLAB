function [R,DER] = expmrot( w )
%
%
%  no chequea nada..!! ni que sea skew symmetric, ni que las dims esten
%  bien, ni que sea cuadrada. nada de nada
%
%
%
%{
if 0

  w = randn(1,1); maxnorm( expmrot( w ) , expm( skewmatrix( w ) ) ); maxnorm( expmrot( skewmatrix(w) ) , expmrot( w ) );
  w = randn(3,1); maxnorm( expmrot( w ) , expm( skewmatrix( w ) ) ); maxnorm( expmrot( skewmatrix(w) ) , expmrot( w ) );
  w = randn(6,1); maxnorm( expmrot( w ) , expm( skewmatrix( w ) ) ); maxnorm( expmrot( skewmatrix(w) ) , expmrot( w ) );

  w = randn(1,1); maxnorm( d_expm( skewmatrix(w) )*d_skewmatrix(w) ,  d_expmrot(w) ); maxnorm( d_expmrot(w) , d_expmrot( skewmatrix(w) )*d_skewmatrix(w) );
  w = randn(3,1); maxnorm( d_expm( skewmatrix(w) )*d_skewmatrix(w) ,  d_expmrot(w) ); maxnorm( d_expmrot(w) , d_expmrot( skewmatrix(w) )*d_skewmatrix(w) );
  w = randn(6,1); maxnorm( d_expm( skewmatrix(w) )*d_skewmatrix(w) ,  d_expmrot(w) ); maxnorm( d_expmrot(w) , d_expmrot( skewmatrix(w) )*d_skewmatrix(w) );


  
  w = randn(1,1); maxnorm( expmrot_mx( w ) , expm( skewmatrix( w ) ) ); maxnorm( expmrot_mx( skewmatrix(w) ) , expmrot_mx( w ) );
  w = randn(3,1); maxnorm( expmrot_mx( w ) , expm( skewmatrix( w ) ) ); maxnorm( expmrot_mx( skewmatrix(w) ) , expmrot_mx( w ) );
  w = randn(6,1); maxnorm( expmrot_mx( w ) , expm( skewmatrix( w ) ) ); maxnorm( expmrot_mx( skewmatrix(w) ) , expmrot_mx( w ) );

  w = randn(1,1); maxnorm( d_expm( skewmatrix(w) )*d_skewmatrix(w) ,  d_expmrot(w) ); maxnorm( d_expmrot(w) , d_expmrot( skewmatrix(w) )*d_skewmatrix(w) );
  w = randn(3,1); maxnorm( d_expm( skewmatrix(w) )*d_skewmatrix(w) ,  d_expmrot(w) ); maxnorm( d_expmrot(w) , d_expmrot( skewmatrix(w) )*d_skewmatrix(w) );
  w = randn(6,1); maxnorm( d_expm( skewmatrix(w) )*d_skewmatrix(w) ,  d_expmrot(w) ); maxnorm( d_expmrot(w) , d_expmrot( skewmatrix(w) )*d_skewmatrix(w) );


  w = randn(6,1); W = skewmatrix( w );
  tic; for it=1:10000, expm(W);       end; toc
  tic; for it=1:10000, expmrot(W);    end; toc
  tic; for it=1:10000, expmrot(w);    end; toc
  tic; for it=1:10000, expmrot_mx(W); end; toc
  tic; for it=1:10000, expmrot_mx(w); end; toc
  




  W = skewmatrix( randn(1,3) );
  d = []; ts = geospace(1e-60,1000,2000);
  for t = ts
    d(end+1)=maxnorm( expm(t*W) - expmrot_mx(t*W) );
  end
  loglog( ts , d , '.' );
  



w = randn(3,1); W = skewmatrix( w );
tic; for it=1:1000, d_expm(W);    end; toc
tic; for it=1:1000, d_expmrot(W); end; toc
tic; for it=1:1000, [a,b]=expmrot(w); end; toc
tic; for it=1:1000, d_expmrot(w); end; toc



  w = normalize( randn(1,3) );
  d = []; ts = geospace(1e-60,1000,2000);
  for t = ts, d(end+1)=maxnorm( d_expm( skewmatrix(t*w) )*d_skewmatrix(w) - d_expmrot(t*w) ); end
  loglog( ts , d , '.' );


end
%}

    asMATRIX = false;
    if ~isvector( w )
        asMATRIX = true;
        switch size(w,1)
            case 2,     w = w(2,1);
            case 3,     w = w([6;7;2]);
        end
    end

    switch numel(w)
        case 1
            c = cos(w);
            s = sin(w);

            R = [ c , -s ; s , c ];

            if nargout > 1
                DER = [ -s ; c ; -c ; -s ];

                if asMATRIX
                    DER = DER * [0,0.5,-0.5,0];
                end

            end

        case 3
            t = sqrt( w(:).' * w(:) );


            if t > 1e-30

                w = w/t;
                x = w(1);
                y = w(2);
                z = w(3);

                s  = sin(t/2);       ss = 2*s*s;
                c  = cos(t/2);       cs = 2*c*s;

                R = [ (x*x-1)*ss + 1  ,  x*y*ss - z*cs   ,  x*z*ss + y*cs   ;...
                    x*y*ss + z*cs   ,  (y*y-1)*ss + 1  ,  y*z*ss - x*cs   ;...
                    x*z*ss - y*cs   ,  y*z*ss + x*cs   ,  (z*z-1)*ss+1    ];

            else

                R = eye(3);

            end

            if nargout > 1

                if t > 1e-30

                    c = cos(t);  s = sin(t);

                    if t > 0.033
                        A  =  - s / t;
                        B  = ( c - 1 )/t;
                    elseif t > 0.003
                        t2 = t*t;
                        A = - 1       + t2/6           - t2*t2/120               + t2*t2*t2/5040;
                        B =     - t/2        + t2*t/24             - t2*t2*t/720                ;
                    else
                        t2 = t*t;
                        A = - 1       + t2/6           - t2*t2/120;
                        B =     - t/2        + t2*t/24            ;
                    end

                    C   = c + A;
                    D   = s + 2 * B;

                    xx  = x*x;
                    xz  = x*z;
                    xy  = x*y;
                    yy  = y*y;
                    yz  = y*z;
                    zz  = z*z;

                    xxx = xx*x;
                    xxy = xx*y;
                    xxz = xx*z;
                    xyy = x*yy;
                    xyz = xy*z;
                    xzz = x*zz;
                    yyy = yy*y;
                    yyz = yy*z;
                    yzz = y*zz;
                    zzz = zz*z;

                    DER = - [   D*x  ,   s*y  ,   s*z  ;...
                        B*y  ,   B*x  ,   A    ;...
                        B*z  , - A    ,   B*x  ;...
                        B*y  ,   B*x  , - A    ;...
                        s*x  ,   D*y  ,   s*z  ;...
                        A    ,   B*z  ,   B*y  ;...
                        B*z  ,   A    ,   B*x  ;...
                        - A    ,   B*z  ,   B*y  ;...
                        s*x  ,   s*y  ,   D*z  ] + ...
                        C*[   0    ,   0    ,   0    ;...
                        xz   ,   yz   ,   zz   ;...
                        - xy   , - yy   , - yz   ;...
                        - xz   , - yz   , - zz   ;...
                        0    ,   0    ,   0    ;...
                        xx   ,   xy   ,   xz   ;...
                        xy   ,   yy   ,   yz   ;...
                        - xx   , - xy   , - xz   ;...
                        0    ,   0    ,   0    ] +...
                        D*[   xxx  ,   xxy  ,   xxz  ;...
                        xxy  ,   xyy  ,   xyz  ;...
                        xxz  ,   xyz  ,   xzz  ;...
                        xxy  ,   xyy  ,   xyz  ;...
                        xyy  ,   yyy  ,   yyz  ;...
                        xyz  ,   yyz  ,   yzz  ;...
                        xxz  ,   xyz  ,   xzz  ;...
                        xyz  ,   yyz  ,   yzz  ;...
                        xzz  ,   yzz  ,   zzz  ];

                else

                    DER = [0 0 0;0,0,1;0,-1,0;0,0,-1;0 0 0;1,0,0;0,1,0;-1,0,0;0 0 0];

                end

                if asMATRIX
                    DER = ( 0.5 * DER ) * [0,0,0,0,0,1,0,-1,0;0,0,-1,0,0,0,1,0,0;0,1,0,-1,0,0,0,0,0];
                end

            end


        otherwise
            if ~asMATRIX
                R = expm( skewmatrix( w ) );
                if nargout > 1
                    DER = d_expm( skewmatrix( w ) )*d_skewmatrix( w );
                end
            else

                R = expm( w );

                if nargout > 1
                    DER = d_expm( w );
                end

            end

    end

end
