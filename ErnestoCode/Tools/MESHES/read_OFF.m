function M = read_OFF( filename )
%
% mesh = read_OFF( filename )
%


%   fid = fopen( filename , 'r' );
%   if( fid == -1 )
%       error('Cant open the file.');
%   end
%   
%   fclose(fid);

  todo = readFile( filename );

  
  i = 1;
  while ~strcmpi( strtrim( todo{i} ) , 'OFF' )
    i = i+1;
    if i > numel( todo )
      error('the doesn''t look like a .OFF file');
    end
  end
  
  i = i+1;
  while isempty( strtrim( todo{i} ) ) || strncmp( todo{i} , '#' , 1 )
    i = i+1;
  end
  
  Nv_Nf = sscanf( todo{i}, '%d' );
  
  M.xyz = zeros( Nv_Nf(1) , 3 );
  for j =  i+1 : i+Nv_Nf(1)
     M.xyz(j-i,:) = sscanf( todo{j}, '%g' );
  end
    
  M.tri = zeros( Nv_Nf(1) , 3 );
  for j = i+Nv_Nf(1)+1 : i+Nv_Nf(1)+ Nv_Nf(2)
    t = sscanf( todo{j}, '%d' );
    M.tri( j - i-Nv_Nf(1) ,: ) = t(2:4);
  end
  M.tri = M.tri + 1;
  
  
end
