function s= SCAM2Mesh( m )
% 
%   s= SCAM2Mesh( mesh )
% 

D= m.SCAM.D;

x= m.SCAM.O(1);
for i=1:m.SCAM.N(1)
  x= [x x(end)+D];
end

y= m.SCAM.O(2);
for i=1:m.SCAM.N(2)
  y= [y y(end)+D];
end

z= m.SCAM.O(3);
for i=1:m.SCAM.N(3)
  z= [z z(end)+D];
end

[x y z]= ndgrid( x,y,z);

s.xyz= [ x(:) y(:) z(:) ];
s.tri= [];

nx= m.SCAM.N(1)+1;
ny= m.SCAM.N(2)+1;
nz= m.SCAM.N(3)+1;

for i=1:nx
  for j=1:ny
    for k=1:nz
      if i<nx
        s.tri(end+1,:)= [ (k-1)*nx*ny+(j-1)*nx+i  (k-1)*nx*ny+(j-1)*nx+i+1 0];
      end
      if j<ny
        s.tri(end+1,:)= [ (k-1)*nx*ny+(j-1)*nx+i  (k-1)*nx*ny+j*nx+i 0];
      end
      if k<nz
        s.tri(end+1,:)= [ (k-1)*nx*ny+(j-1)*nx+i  k*nx*ny+(j-1)*nx+i 0];
      end
    end
  end
end




