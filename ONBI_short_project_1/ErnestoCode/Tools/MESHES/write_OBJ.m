function write_OBJ( mesh , filename )
%
%  write_OBJ( 'filename' , mesh );
%

fid = fopen(filename,'wt');
if( fid==-1 )
    error('Cant open the file.');
end

mesh= FixMesh( mesh );

% header
fprintf(fid, '# number of points : %d \n', size(mesh.xyz,1) );
fprintf(fid, '# number of faces  : %d \n', size(mesh.tri,1));

% write the points & faces
fprintf(fid, '\n# points\n');
fprintf(fid, 'v %.10e %.10e %.10e\n', mesh.xyz');

if isfield( mesh, 'uv')
  fprintf(fid, '\n# texture map\n');
  fprintf(fid, 'vt %.10e %.10e\n', mesh.uv');
end

if isfield( mesh, 'xyzNORMALS')
  fprintf(fid, '\n# normals \n');
  fprintf(fid, 'vn %.10e %.10e %.10e\n', mesh.xyzNORMALS');
end


fprintf(fid, '\n# faces\n');
if      (  isfield( mesh, 'uv') &&  isfield( mesh, 'xyzNORMALS') )
  fprintf(fid, 'f %d/%d/%d %d/%d/%d %d/%d/%d\n', mesh.tri(:,[1 1 1 2 2 2 3 3 3])');
elseif  ( ~isfield( mesh, 'uv') &&  isfield( mesh, 'xyzNORMALS') )
  fprintf(fid, 'f %d//%d %d//%d %d//%d\n', mesh.tri(:,[1 1 2 2 3 3])');
elseif  (  isfield( mesh, 'uv') && ~isfield( mesh, 'xyzNORMALS') )
  fprintf(fid, 'f %d/%d %d/%d %d/%d\n', mesh.tri(:,[1 1 2 2 3 3])');
else
  fprintf(fid, 'f %d %d %d\n', mesh.tri');
end

fclose(fid);
