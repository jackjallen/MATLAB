function m= MeshContour( D, X, Y, Z , value)
% 
% m= MeshContour( D, X, Y, Z , value)
% 

if nargin < 5
  value = 0;
end

[f,v]= isosurface( X,Y,Z,D, value );
m.xyz= v;
m.tri= f;

