function mesh= CreateLocators( mesh , n )
%
% mesh= CreateLocators( mesh , n )
%
% create a structured cartesian auxiliar mesh
% to be used by the locators and fill the
% nodes and elements locators structures.
%
% n is the number of buckets in the
%    SCAM structure
%
% if n is not specified, them, automatic size adopted
%
% if n is equal to zero the auxiliar structures are
%    removed.
%

if nargin < 2, n= -1; end

if n ~= 0
  minx= min( mesh.xyz(:,1) );
  maxx= max( mesh.xyz(:,1) );
  miny= min( mesh.xyz(:,2) );
  maxy= max( mesh.xyz(:,2) );
  minz= min( mesh.xyz(:,3) );
  maxz= max( mesh.xyz(:,3) );

  dx= maxx-minx;
  dy= maxy-miny;
  dz= maxz-minz;

  if n < 0
    a= sqrt( (dx*dy+dx*dz+dy*dz)/size( mesh.tri , 1 ) );
    n= ceil(max([dx dy dz])/a/10);
  end

  d= max( [dx dy dz] )/n;

  center= [ minx+maxx , maxy+miny , maxz+minz ]/2;

  nx= ceil( dx/d ); if ~nx, nx=1; end;
  ny= ceil( dy/d ); if ~ny, ny=1; end;
  nz= ceil( dz/d ); if ~nz, nz=1; end;

  mesh.SCAM.D = d;
  mesh.SCAM.O = center - d/2*[nx ny nz];
  mesh.SCAM.N = [nx ny nz];

  mesh= CreateELOC( CreatePLOC( mesh ));
else
  if isfield(mesh,'SCAM')
    mesh=rmfield(mesh,'SCAM');
  end
  if isfield(mesh,'ELOC')
    mesh=rmfield(mesh,'ELOC');
  end
  if isfield(mesh,'PLOC')
    mesh=rmfield(mesh,'PLOC');
  end
end

