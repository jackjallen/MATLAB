function write_VTP( M , filename , ab )
%
% write_VTP( mesh , filename , ['ascii' or 'binary'])
%

  if nargin < 3, ab = 'ascii'; end
  ab = upper(ab);
  if ~strcmp( ab , 'ASCII' ) && ~strcmp( ab , 'BINARY' )
    error( '''ascii'' or ''binary'' were expected');
  end

  fid = fopen(filename,'w');
  if( fid==-1 ), error('Cant open the file.'); end
  CLEANUP = onCleanup( @()fclose(fid) );
  

  % header
  fprintf(fid, '# vtk DataFile Version 3.0\n');
  fprintf(fid, 'vtk output\n');
  switch ab
    case 'ASCII',  fprintf(fid, 'ASCII\n');
    case 'BINARY', fprintf(fid, 'BINARY\n');
  end
  fprintf(fid, 'DATASET POLYDATA\n');

  % points
  N_xyz = 0;
  if isfield(M,'xyz') && size( M.xyz , 1 ) > 0
    N_xyz = size( M.xyz , 1 );
    switch ab
      case 'ASCII',   fprintf(fid, '\nPOINTS %d double\n', N_xyz );
      case 'BINARY',  fprintf(fid, '\nPOINTS %d %s\n', N_xyz , class_as_c( M.xyz ) );
    end
    write_in_file(fid, M.xyz );
  end


  % polygons
  N_tri = 0;
  if isfield( M , 'tri' ) && size( M.tri , 1 ) > 0
    M.tri = int32( M.tri.' - 1 );
    M.tri = M.tri( : , ~all( M.tri < 0 , 1 ) );

    N_tri = size( M.tri , 2 );
    
    switch ab
      case 'ASCII'
        
        M.tri = [ int32( sum( M.tri >= 0 , 1 ) ) ; M.tri ; ones( 1 , N_tri , 'int32' )*intmin('int32') ];
        M.tri = M.tri(:);
        M.tri( M.tri < 0 & M.tri > intmin('int32') ) = [];
        M.tri( M.tri == intmin('int32') ) = -1;

        str = sprintf('%d ', M.tri );
        str = strrep( str , ' -1 ' , char(10) );
        str(end) = [];
        fprintf( fid , 'POLYGONS %d %d\n', N_tri , numel( M.tri ) - N_tri );
        fprintf( fid , str );
        
      case 'BINARY'

        M.tri = [ int32( sum( M.tri >= 0 , 1 ) ) ; M.tri ];
        M.tri = M.tri(:);
        M.tri( M.tri < 0 ) = [];

        fprintf( fid , 'POLYGONS %d %d\n', N_tri , numel( M.tri ) );
%         fwrite_flipped( fid , M.tri );
        fwrite( fid , M.tri , 'int32' , 0 , 'b' );

    end
    
    fprintf(fid, '\n');
  end



  fields = setdiff( fieldnames( M ) , {'xyz','tri','xyzNORMALS','triNORMALS','uv'} );

  % points data
  if N_xyz
    fprintf( fid , '\nPOINT_DATA %d\n', N_xyz );
    
    if isfield( M , 'uv' )
      switch ab
        case 'ASCII',  fprintf( fid , '\nTEXTURE_COORDINATES UV 2 double\n');
        case 'BINARY', fprintf( fid , '\nTEXTURE_COORDINATES UV 2 %s\n', class_as_c( M.uv ) );
      end
      write_in_file(fid, M.uv );
    end

    if isfield( M , 'xyzNORMALS')
      switch ab
        case 'ASCII',  fprintf( fid , '\nNORMALS Normals double\n');
        case 'BINARY', fprintf( fid , '\nNORMALS Normals %s\n', class_as_c( M.xyzNORMALS ) );
      end
      write_in_file(fid, M.xyzNORMALS );
    end

    for f = 1:numel(fields)
      field= fields{f};
      if strncmp( field, 'xyz',3)
        fprintf(fid, '\nFIELD field 1\n');
        switch ab
          case 'ASCII',  fprintf(fid, '%s %d %d double\n', field(4:end), size(M.(field),2), N_xyz );
          case 'BINARY', fprintf(fid, '%s %d %d %s\n', field(4:end), size(M.(field),2), N_xyz , class_as_c( M.(field) ) );
        end
        write_in_file( fid, M.(field) );
      end
    end
  end

  % cells data
  if N_tri
    fprintf(fid,'\nCELL_DATA %d', N_tri );

    if isfield( M , 'triNORMALS')
      switch ab
        case 'ASCII',  fprintf( fid , '\nNORMALS CellsNormals double\n');
        case 'BINARY', fprintf( fid , '\nNORMALS CellsNormals %s\n', class_as_c( M.triNORMALS ) );
      end
      write_in_file(fid, M.triNORMALS );
    end

    for f = 1:numel(fields)
      field= fields{f};
      if strncmp( field, 'tri',3)
        fprintf(fid, '\nFIELD field 1\n');
        switch ab
          case 'ASCII',  fprintf(fid, '%s %d %d double\n', field(4:end), size(M.(field),2), N_tri );
          case 'BINARY', fprintf(fid, '%s %d %d %s\n', field(4:end), size(M.(field),2), N_tri , class_as_c( M.(field) ) );
        end
        write_in_file( fid, M.(field) );
      end
    end
  end

  fprintf( fid, '\n' );


  function c = class_as_c( x )
    switch class(x)
      case 'double',  c = 'double';
      case 'single',  c = 'float';
      case 'uint64',  c = 'unsigned_long';
      case 'int64',   c = 'long';
      case 'uint32',  c = 'unsigned_int';
      case 'int32',   c = 'int';
      case 'uint16',  c = 'unsigned_short';
      case 'int16',   c = 'short';
      case 'uint8',   c = 'unsigned_char';
      case 'int8',    c = 'char';
      case 'logical', c = 'bool';
    end
  end

  function o = fwrite_flipped( fid , x )
    x1 = x(1);
    x = x(:);
    x = typecast( x , 'uint8' );
    x = reshape( x , numel(typecast(x1,'uint8')) , [] );
    x = flipdim( x , 1 );
    x = x(:);
    o = fwrite( fid , x , 'uint8' );
  end
  function str = mat_2_str( x )

%         fprintf(fid, '%.10e %.10e %.10e\n', M.xyz.' );

        %{
        str = regexprep(regexprep(regexprep(regexprep(regexprep( ...
                  sprintf('%.25e %.25e %.25e\n',M.xyz.' ) ...
              ,'e+000',''),'0*e','e'),'\.e','e'),'e\+0*','e'),'e-0*','e-');
        %}
        %{
        str = strrep( ...
                  sprintf(' %0.025f %0.025f %0.025f \n',M.xyz.' ) ...
              ,'0 ', ' ' , ' 0' ,' ');
        %}
%         str = regexprep( regexprep( sprintf(' %0.025f %0.025f %0.025f \n', M.xyz.' ) ,'0* ', ' ') , ' 0*' , ' ');
        
%         str = sprintf(' %0.025g %0.025g %0.025g \n', [10 10000 0.00001; 0 pi 1e-38].' );
%         str = sprintf(' %0.025e %0.025e %0.025e \n', [ [10 10000 0.00001; 0 pi 1e-38].' , single(M.xyz(1:20,:).') ] );
    
    
    str = sprintf( [ repmat( ' % 0.018e  ' , 1 , size( x , 2 ) ) , ' \n' ] , x' );

%     str_n = -1; while str_n ~= numel(str)
%       str_n = numel(str);
%       str = builtin( 'strrep' , str , '  ' , ' ' );
%     end
%     str = builtin( 'strrep' , str , 'e+000' , ' '  );
%     str = builtin( 'strrep' , str , 'e+00'  , 'e+' );
%     str = builtin( 'strrep' , str , 'e+0'   , 'e+' );
%     str = builtin( 'strrep' , str , 'e+'    , 'e'  );
%     str = builtin( 'strrep' , str , 'e-00'  , 'e-' );
%     str = builtin( 'strrep' , str , 'e-0'   , 'e-' );
%     str_n = -1; while str_n ~= numel(str)
%       str_n = numel(str);
%       str = builtin( 'strrep' , str , '0e' , 'e' );
%     end
%     str = builtin( 'strrep' , str , '.e'   , 'e' );
%     str_n = -1; while str_n ~= numel(str)
%       str_n = numel(str);
%       str = builtin( 'strrep' , str , '00 ' , '0 ' );
%     end
%     str_n = -1; while str_n ~= numel(str)
%       str_n = numel(str);
%       str = builtin( 'strrep' , str , ' 00' , ' 0' );
%     end
%     str = builtin( 'strrep' , str , ' 0.0 ' , ' 0 ' );
%     str = builtin( 'strrep' , str , ' .' , ' 0.' );
%     str = builtin( 'strrep' , str , '. ' , ' ' );

    str = builtin( 'strrep' , str , [' ' char(10)] , char(10) );
    str = builtin( 'strrep' , str , [char(10) ' '] , char(10) );
    if str(1) == ' ', str = str(2:end); end
    
  end

  function write_in_file( fid , x )
    switch ab
      case 'ASCII'
%         fprintf( fid , mat_2_str( x ) );
        str = sprintf( [ repmat( ' % 0.018e  ' , 1 , size( x , 2 ) ) , ' \n' ] , x' );
        fprintf( fid , str , x.' );
        
      case 'BINARY'
%         fwrite_flipped( fid , x.' );
        fwrite( fid , x.' , class(x) , 0 , 'b' );
    end
    fprintf(fid, '\n');
  end


end
