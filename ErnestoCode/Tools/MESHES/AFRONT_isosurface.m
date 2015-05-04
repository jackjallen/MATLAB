function M = AFRONT_isosurface( X , value , varargin )
% 
% A wrapper to afront isosurface
% 

  if ispc
    exe = fullfile( fileparts( which( mfilename )) , 'AFRONT' , 'afront.exe' );
  elseif isunix
    exe = fullfile( fileparts( which( mfilename )) , 'AFRONT' , 'afront_0_2' );
  end

  CLEANUP = {};
  
  dirname = tmpname('AFRONT_isosurface_');
  mkdir( dirname );
  CLEANUP = [ onCleanup( @() rmdir(dirname,'s') ) , CLEANUP ];
  
  cdw = pwd;
  CLEANUP = [ onCleanup( @() cd(cdw) ) , CLEANUP ];
  
  cd( dirname );
  
  
  
  in_data  = fullfile(  dirname , 'data.nhdr' );
  out_mesh = fullfile(  dirname , 'mesh.m'    );

  SPATIAL_TRANSFORM = [];
  if isI3D( X )
    X.data = single( X.data );
    X = crop( X , 2 ,  imoutline( X <= value/1.1 ) | imoutline( X >= value*1.1 ) );
    
    X = X.matrix2coords;
    X = rot90( X , rot90( X , eye(3) ) );
    SPATIAL_TRANSFORM = X.SpatialTransform;
    X.SpatialTransform = eye(4);
  end
  
  write_NRRD( X , in_data );

  command = sprintf('%s  -nogui  %s -outname %s' , exe , in_data , out_mesh );
  
  for i = 1:numel(varargin)
    command = [ command  '  '  varargin{i} ];
  end

%   command = [ command , ' ' , '-marchingcubes ' , sprintf('%.15g',value) ' ' ];
  command = [ command , ' ' , '-tri ' , sprintf('%.15g',value) ' ' ];
%   command = [ command , ' ' , '-tri ' , sprintf('%.15g',value) ' bspline' ];
  
  fprintf('Running Afront ... wait\n' );

%   system( command );
  [o1,o2] = system( command );
%   if o1 ~= 0, disp( o2 ); end
  
  fprintf('OK\n');
  
  M = read_MESH( out_mesh );
  
  if ~isempty( SPATIAL_TRANSFORM )
    M.xyz = transform( M.xyz , SPATIAL_TRANSFORM );
  end

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