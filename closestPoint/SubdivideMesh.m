function m= SubdivideMesh( m , repeat )
% 
% m= Subdivide( m , times )
% 

  if nargin < 2
    repeat = 1;
  end

  for k=1:repeat
    n_tri= size( m.tri,1 );
    n_xyz= size( m.xyz,1 );

    m.xyz( n_xyz+3*n_tri,3)=0;
    m.tri( n_tri+3*n_tri,3)=1;

    idp= n_xyz;
    idc= n_tri;
    for t=1:n_tri
      c= m.tri(t,:);
      p1= m.xyz( c(1) , :);
      p2= m.xyz( c(2) , :);
      p3= m.xyz( c(3) , :);
      m.xyz( idp+1,: )= (p1+p2)/2;
      m.xyz( idp+2,: )= (p2+p3)/2;
      m.xyz( idp+3,: )= (p3+p1)/2;
      m.tri( t,: )= [ c(1) idp+1 idp+3 ]; 
      m.tri( idc+1,: )=[ c(2) idp+2 idp+1 ];
      m.tri( idc+2,: )=[ c(3) idp+3 idp+2 ];
      m.tri( idc+3,: )=[ idp+1 idp+2 idp+3 ];
      idp= idp+3;
      idc= idc+3;
    end
    
    try, m = TidyMesh( m ); end
  end

%   m= CleanMesh( m );

end
