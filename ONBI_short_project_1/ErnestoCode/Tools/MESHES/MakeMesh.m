function mesh= MakeMesh( varargin )
% 
% mesh= MakeMesh( 'Points', xyz
%                 'Triangles', tri
%
         

zz= parseargs( varargin,'points','p');
if zz
  mesh.xyz= varargin{zz+1};
end
if size( mesh.xyz,2 ) < 3
  mesh.xyz(:,3)=0;
end

zz= parseargs( varargin,'triangles','tri','t','triangle');
if zz
  mesh.tri= varargin{zz+1};
else
  mesh.tri= delaunayn( mesh.xyz(:,[1 2 3]) );
end
