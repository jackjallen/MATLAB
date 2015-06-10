function M = read_VTK_POLYDATA( filename )
%
%


  fid = fopen(filename,'r');
  if( fid == -1 ), error('Cant open the file.'); end
  CLEANUP = onCleanup( @()fclose(fid) );
  
  ab = 'ASCII';
  
  nowREADING = '';
  while ~feof(fid)
    l = fgetl( fid );
    [key,l] = strtok( l );
    key = upper( key );
    switch key
        case '#',       continue;
        case 'BINARY',  ab = 'BINARY';
        case 'ASCII',   ab = 'ASCII';
        case 'POLYDATA'
            dataset = strtok( l );
            if ~strcmp( dataset , 'POLYDATA' )
                error('esto esta preparado para leer POLYDATA solamente');
            end
            continue;
        case 'POINTS'
            points = textscan( l , '%f %s');
            READnum  = points{1}*3;
            READtype = class_as_matlab( points{2}{1} );
            M.xyz = read(); 
            M.xyz = reshape( M.xyz , [ 3 , numel(M.xyz)/3 ] ).';
            M.xyz = double(M.xyz);

        case 'POLYGONS'
            READtype = 'int32';
            cells = textscan( l , '%f %f');
            READnum = cells{2};
            CELLSnum = cells{1};
            M.tri = read();
            
            M.tri = reshape( M.tri , [ READnum/CELLSnum , CELLSnum ] ).';
            M.tri(:,1) = [];
            M.tri = M.tri + 1;
            M.tri = double(M.tri);
            break;
    end
      
  end
  
  
    function x = read()
       switch ab
           case 'BINARY'
               x = fread( fid , READnum , ['*' READtype] );
               
           case 'ASCII'
               x = textscan( fid , '%f' , READnum );
               x = x{1};
               
       end
    end
  
  
  
%   % polygons
%   N_tri = 0;
%   if isfield( M , 'tri' ) && size( M.tri , 1 ) > 0
%     M.tri = int32( M.tri.' - 1 );
%     M.tri = M.tri( : , ~all( M.tri < 0 , 1 ) );
% 
%     N_tri = size( M.tri , 2 );
%     
%     switch ab
%       case 'ASCII'
%         
%         M.tri = [ int32( sum( M.tri >= 0 , 1 ) ) ; M.tri ; ones( 1 , N_tri , 'int32' )*intmin('int32') ];
%         M.tri = M.tri(:);
%         M.tri( M.tri < 0 & M.tri > intmin('int32') ) = [];
%         M.tri( M.tri == intmin('int32') ) = -1;
% 
%         str = sprintf('%d ', M.tri );
%         str = strrep( str , ' -1 ' , char(10) );
%         str(end) = [];
%         fprintf( fid , '\nCELLS %d %d\n', N_tri , numel( M.tri ) - N_tri );
%         fprintf( fid , str );
%         
%       case 'BINARY'
% 
%         M.tri = [ int32( sum( M.tri >= 0 , 1 ) ) ; M.tri ];
%         M.tri = M.tri(:);
%         M.tri( M.tri < 0 ) = [];
% 
%         fprintf( fid , '\nCELLS %d %d\n', N_tri , numel( M.tri ) );
% %         fwrite_flipped( fid , M.tri );
%         fwrite( fid , M.tri , 'int32' , 0 , 'b' );
% 
%     end
%     
%     fprintf(fid, '\n');
%   end
% 
%   %CELL_TYPE
%   if isscalar( M.celltype )
%       M.celltype = zeros( [ N_tri , 1 ] , 'int32' ) + int32( M.celltype );
%   elseif numel( M.celltype ) ~= N_tri
%       error('numel of celltype does not coincide with number of cells');
%   else
%       M.celltype = int32( M.celltype );
%   end
%   fprintf(fid,'\nCELL_TYPES %d\n',N_tri);
%   write_in_file(fid, M.celltype );
%   
%   
%   
% 
%   fields = setdiff( fieldnames( M ) , {'xyz','tri','xyzNORMALS','triNORMALS','uv'} );
% 
%   % points data
%   if N_xyz
%     fprintf( fid , '\nPOINT_DATA %d\n', N_xyz );
%     
%     if isfield( M , 'uv' )
%       switch ab
%         case 'ASCII',  fprintf( fid , '\nTEXTURE_COORDINATES UV 2 double\n');
%         case 'BINARY', fprintf( fid , '\nTEXTURE_COORDINATES UV 2 %s\n', class_as_c( M.uv ) );
%       end
%       write_in_file(fid, M.uv );
%     end
% 
%     if isfield( M , 'xyzNORMALS')
%       switch ab
%         case 'ASCII',  fprintf( fid , '\nNORMALS Normals double\n');
%         case 'BINARY', fprintf( fid , '\nNORMALS Normals %s\n', class_as_c( M.xyzNORMALS ) );
%       end
%       write_in_file(fid, M.xyzNORMALS );
%     end
% 
%     for f = 1:numel(fields)
%       field= fields{f};
%       if strncmp( field, 'xyz',3)
%         fprintf(fid, '\nFIELD field 1\n');
%         switch ab
%           case 'ASCII',  fprintf(fid, '%s %d %d double\n', field(4:end), size(M.(field),2), N_xyz );
%           case 'BINARY', fprintf(fid, '%s %d %d %s\n', field(4:end), size(M.(field),2), N_xyz , class_as_c( M.(field) ) );
%         end
%         write_in_file( fid, M.(field) );
%       end
%     end
%   end
% 
%   % cells data
%   if N_tri
%     fprintf(fid,'\nCELL_DATA %d', N_tri );
% 
%     if isfield( M , 'triNORMALS')
%       switch ab
%         case 'ASCII',  fprintf( fid , '\nNORMALS CellsNormals double\n');
%         case 'BINARY', fprintf( fid , '\nNORMALS CellsNormals %s\n', class_as_c( M.triNORMALS ) );
%       end
%       write_in_file(fid, M.triNORMALS );
%     end
% 
%     for f = 1:numel(fields)
%       field= fields{f};
%       if strncmp( field, 'tri',3)
%         fprintf(fid, '\nFIELD field 1\n');
%         switch ab
%           case 'ASCII',  fprintf(fid, '%s %d %d double\n', field(4:end), size(M.(field),2), N_tri );
%           case 'BINARY', fprintf(fid, '%s %d %d %s\n', field(4:end), size(M.(field),2), N_tri , class_as_c( M.(field) ) );
%         end
%         write_in_file( fid, M.(field) );
%       end
%     end
%   end
% 
%   fprintf( fid, '\n' );


  function c = class_as_matlab( x )
    switch lower(strtrim(x))
      case 'double',         c = 'double';
      case 'float',          c = 'single';
      case 'unsigned_long',  c = 'uint64';
      case 'long',           c = 'int64';
      case 'unsigned_int',   c = 'uint32';
      case 'int',            c = 'int32';
      case 'unsigned_short', c = 'uint16';
      case 'short',          c = 'int16';
      case 'unsigned_char',  c = 'uint8';
      case 'char',           c = 'int8';
    end
  end



end