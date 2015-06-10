function M = AFRONT_Remesh( M , varargin )
% 
% A wrapper to afront remeshing
% 

  if numel(varargin) > 0  &&  isstruct( varargin{1} )
    M2 = varargin{1};
    varargin(1) = [];
  else
    M2 = [];
  end
  M2 = [];  %% descomentar esto si entiendo como funciona csg
    
  location = fullfile( fileparts( which('AFRONT_Remesh.m') ) , 'AFRONT' , 'afront.exe' );
  
  
  in_mesh  = [ tempname '.off' ];
  in_mesh2 = [ tempname '.off' ];

  out_mesh = [ tempname '.m'   ];
  
  write_OFF( M , in_mesh );
  if ~isempty( M2 )
    write_OFF( M2 , in_mesh2 );
  end
  
  
  if isempty( M2 )
    command = sprintf('%s  -nogui  %s -outname %s' , location , in_mesh , out_mesh );
  else
    command = sprintf('%s  -nogui  %s %s -outname %s' , location , in_mesh , in_mesh2 , out_mesh );
  end
  
  for i = 1:numel(varargin)
    command = [ command  '  '  varargin{i} ];
  end
  
  command = [ command '  ' '-tri' ];
  
  fprintf('Running Afront ... wait\n' );

  system( command );
%   [o1,o2] = system( command );
%   if o1 ~= 0, disp( o2 ); end
  
  fprintf('OK\n');
  
  M = read_MESH( out_mesh );

end


%{

cd d:\afront\
fid = fopen( 'iso.nhdr.raw' , 'r' )
xx = fread( fid , 'uchar=>double' );
xx = reshape( xx, [256 256 256]);

M = FixMesh( isosurface( xx(90:140,65:140,140:240) , 200 ) )
M = vtkPolyDataConnectivityFilter( M )
M, M = CleanMesh( M , 1e-4 )
plotMESH( M )

Mr = AFRONT_Remesh( M ,'-rho .5 -min_step .5' )
plotMESH( Mr )

%}