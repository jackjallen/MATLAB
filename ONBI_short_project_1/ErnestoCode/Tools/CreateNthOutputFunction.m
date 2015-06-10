function  CreateNthOutputFunction( fun , N )

  if ~ischar( fun ), error('an string expected'); end
  if nargin < 2 , N = 10; end

  filetype = exist( fun );
  if ~any( filetype == [2 3 5 6] ), error('invalid function'); end
  
  try
    nouts = nargout( fun );
  catch
    nouts = -1;
  end
  if nouts == -1, nouts = N; end


  p = fullfile( fileparts( mfilename('fullpath') ) , 'NthOutputFunctions' );
  if ~isdir( p ), mkdir( p ); end


  [kk,funname,kk] = fileparts( fun );
  
  for no = 1:nouts
    
    fname = fullfile( p , funname );
    fname = sprintf('%s__o%d.m', fname , no );
    
    fid = fopen( fname , 'w');
    fprintf( fid , 'function out = %s__o%d( varargin )\n' , funname , no );
    fprintf( fid , '  [' );
    for i = 1:no-1    , fprintf(fid,' kk ,'); end
    fprintf( fid , ' out ' );
    for i = no+1:nouts, fprintf(fid,', kk '); end
    fprintf( fid , '] = ' );
    fprintf( fid , '%s( varargin{:} );\n\n' , funname );
    
    fclose( fid );
  end
    
    
  addpath( p );
  
end
