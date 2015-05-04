function M = expmsim( M )
%{ 
tau = @(p) [ p(4) -p(1) p(2) p(5) ; p(1) p(4) -p(3) p(6) ; -p(2) p(3) p(4) p(7) ; 0 0 0 0 ];
t = tau( rand(7,1) );

tic; for i=1:10000, eM  = expm( t ); end, toc
  
tic; for i=1:10000, eMs = expmsim( t ); end, toc

eM - eMs
%} 

  if  ndims( M ) ~= 2  ||  size(M,1) ~= 4  ||  size(M,2) ~= 4 || ...
      M(1) ~= M(6)     ||  M(1) ~= M(11)   ||  ...
      M(2) ~= -M(5)    ||  M(3) ~= -M(9)   ||  M(7) ~= -M(10)  || ...
      any( M(4,:) ~= 0 )
    warning('EXPMSIM:invalidMatrix','The matrix is not tangent to the similarity group. Calling expm .');
    M = expm( M );
    return;
  end

  e   = M(1);
  
  rz  = M(2);
  rx  = M(7);
  ry  = M(9);

  rxx = rx*rx;
  ryy = ry*ry;
  rzz = rz*rz;
  
  rxy = rx*ry;
  rxz = rx*rz;
  ryz = ry*rz;
  
  
  aa  = rxx + ryy + rzz;
  a   = sqrt( aa );

  tx = M(13);
  ty = M(14);
  tz = M(15);
  
  rt  = rx*tx + ry*ty + rz*tz;

  expe = exp(e);
  cosa = cos(a);
  sina = sin(a);

  if a > 1e-7
    cc   = (cosa-1)/aa;
    ss   = sina/a;
  else
    cc   = -1/2 + aa/24;
    ss   =  1 - aa/6;
  end  
  
  
  ee = e*e;
  
  if abs(e) < 1e-8   &&  abs(a) < 1e-8  
    C1 = (2 + e)/2;
    C2 = (-3 - 2*e)/6;
    C3 = 1/6 + e/8;
  elseif abs(e) < 3e-4   &&  abs(a) < 1e-3
    C1 = (20*(6 + e*(3 + e)) - aa*(20 + 3*e*(5 + 2*e)))/120;
    C2 = -1/2 - e/3 + (aa*(15 + e*(12 + 5*e)))/360 - ee/8;
    C3 = 1/6 + (e*(5 + 2*e))/40 - (aa*(42 + 5*e*(7 + 3*e)))/5040;
  elseif abs(e) < 3e-4   &&  abs(a) >= 1e-3
    C1 = (2*a*e*(-1 + cosa*(1 + e)) + (aa*(2 + e*(2 + e)) - 2*ee)*sina)/(2*a^3);
    C2 = (cosa*(aa*(2 + e*(2 + e)) - 2*ee) - 2*(aa - ee + a*e*(1 + e)*sina))/(2*a^4);
    C3 = (6*a*e - 6*a*cosa*e*(1 + e) + a^3*(6 + e*(3 + e)) - 3*(aa*(2 + e*(2 + e)) - 2*ee)*sina)/(6*a^5);
    
%     C1 = ((-1 + cosa)*e + a*(1 + e)*sina)/a^2;
%     C2 = -((a - a*cosa*(1 + e) + e*sina)/a^3);
%     C3 = (2*e - 2*cosa*e + aa*(2 + e) - 2*a*(1 + e)*sina)/(2*a^4);
    
    
  elseif abs(e) >= 3e-4   &&  abs(a) < 1e-3

    C1 = (2*ee*(-1 + expe) + aa*(2 - (2 + (-2 + e)*e)*expe))/(2*ee*e);
    C2 = (-6*ee*(1 + (-1 + e)*expe) + aa*(6 + (-6 + e*(6 + (-3 + e)*e))*expe))/(6*ee*ee);
    C3 = (12*ee*(-2 + (2 + (-2 + e)*e)*expe) + aa*(24 - (24 + e*(-24 + e*(12 + (-4 + e)*e)))*expe))/(24*ee*ee*e);
    
%     C1 = (-1 + expe)/e;
%     C2 = -((1 + (-1 + e)*expe)/ee);
%     C3 = (-2 + (2 + (-2 + e)*e)*expe)/(2*ee*e);
    
  else

    iaaee = 1/(aa+ee);
    C1 = (expe*(cosa*e + a*sina) - e)*iaaee;
    C2 = (expe*(a*cosa - e*sina) - a)*iaaee/a;
    C3 = (expe*(aa + ee - e*(cosa*e + a*sina)) - aa)*iaaee/(aa*e);
    
  end
 
  
  M = [  expe*( cosa  - cc*rxx ) ,  -expe*( ss*rz + cc*rxy )  ,  expe*( ss*ry - cc*rxz )  ,  C1*tx + C2*( rz*ty - ry*tz ) + C3*rt*rx ;... 
         expe*( ss*rz - cc*rxy ) ,   expe*( cosa  - cc*ryy )  , -expe*( ss*rx + cc*ryz )  ,  C1*ty + C2*( rx*tz - rz*tx ) + C3*rt*ry ;...
        -expe*( ss*ry + cc*rxz ) ,   expe*( ss*rx - cc*ryz )  ,  expe*( cosa  - cc*rzz )  ,  C1*tz + C2*( ry*tx - rx*ty ) + C3*rt*rz  ;...
        0 , 0 , 0 , 1];
  

end
