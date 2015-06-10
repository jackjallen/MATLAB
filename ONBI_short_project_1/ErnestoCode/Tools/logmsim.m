function L = logmsim( M )

  s = det( M );
  if s < 0
    error('M has to be a right transform');
  end
  
  
 L = logm( M );
return;

%{
tita = acos((trace(R)-1)/2);
(R-R.')/(2*sin(tita))*tita
%}
 
  s = s^(1/3);
  
  R = M(1:3,1:3) / s;
  
  s = log( s );
  

  if sum(sum(abs( R*R' - eye(3)))) > 1e-10, error('ERROR: R is not a rotation matrix.' ); end
  if abs(det(R) - 1)>1e-10, error('ERROR: R is not a rotation matrix.' ); end
  
  tr = 1 - R(1) - R(5) - R(9);

  sg = sign( R(6) - R(8) );  if sg==0, sg=1; end
  q(1) = abs( sqrt( 2*R(1) + tr )/2 )*sg;

  sg = sign( R(7) - R(3) );  if sg==0, sg=1; end
  q(2) = abs( sqrt( 2*R(5) + tr )/2 )*sg;

  sg = sign( R(2) - R(4) );  if sg==0, sg=1; end
  q(3) = abs( sqrt( 2*R(9) + tr )/2 )*sg;


  nq  = norm( q );
  if nq > 0
    xyz = 2*q/nq*asin( nq );
  
    L = [ s -xyz(3) xyz(2) 0; xyz(3) s -xyz(1) 0; -xyz(2) xyz(1) s 0; 0 0 0 0 ];
    L(1:3,4) = L(1:3,1:3) * ( ( M(1:3,1:3) - eye(3) ) \ M(1:3,4) );
  else
    L = [ s 0 0 M(1,4) ; 0 s 0 M(2,4) ; 0 0 s M(3,4) ; 0 0 0 0 ];
  end



end
