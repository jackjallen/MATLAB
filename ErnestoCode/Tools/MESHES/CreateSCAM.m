function mesh= CreateSCAM( mesh , nx , ny , nz )
%  
% mesh= CreateSCAM( mesh , nx , ny , nz )
% 
% create a structured cartesian auxiliar mesh 
% to be used by the locators.
% 
% If the data XYZ in mesh are 2D (in any dimension!)
% this dimension will be ignored
% 

if nargin < 3
  ny= nx;
  nz= nx;
end

mesh.SCAM.Ox= min( mesh.xyz(:,1) );
mesh.SCAM.Oy= min( mesh.xyz(:,2) );
mesh.SCAM.Oz= min( mesh.xyz(:,3) );

Mx= max( mesh.xyz(:,1) );
My= max( mesh.xyz(:,2) );
Mz= max( mesh.xyz(:,3) );

if ( mesh.SCAM.Ox == Mx )
  mesh.SCAM.Nx = 1;
  mesh.SCAM.Dx = 1;
else
  mesh.SCAM.Nx = nx;
  mesh.SCAM.Dx = ( Mx - mesh.SCAM.Ox )/nx;
end

if ( mesh.SCAM.Oy == My )
  mesh.SCAM.Ny = 1;
  mesh.SCAM.Dy = 1;
else
  mesh.SCAM.Ny = ny;
  mesh.SCAM.Dy = ( My - mesh.SCAM.Oy )/ny;
end

if ( mesh.SCAM.Oz == Mz )
  mesh.SCAM.Nz = 1;
  mesh.SCAM.Dz = 1;
else
  mesh.SCAM.Nz = nz;
  mesh.SCAM.Dz = ( Mz - mesh.SCAM.Oz )/nz;
end


mesh.SCAM.D = [ mesh.SCAM.Dx mesh.SCAM.Dy mesh.SCAM.Dz];
mesh.SCAM.N = [ mesh.SCAM.Nx mesh.SCAM.Ny mesh.SCAM.Nz];
mesh.SCAM.O = [ mesh.SCAM.Ox mesh.SCAM.Oy mesh.SCAM.Oz];
