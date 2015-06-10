function dups = findDuplicates( dir )
  if nargin<1, dir = pwd; end

  files = rdir(dir);
  for i=1:numel(files), dups{i} = { files(i).name }; end
  
  
  fprintf('%6d  of  %6d' , 0 , numel(files) );
  for i= 1:numel(files)
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b');
    fprintf('%6d  of  %6d',i,numel(files));
    
    F = fopen( files(i).name , 'rb' );
    if F > 0
      try
        Fcontent = fread(F , 'int8=>int8' );
        fclose( F );
      catch
        fclose( F );
        continue;
      end
    else
      continue;
    end
    for n= i+1:numel(files)
      if files(i).bytes ~= files(n).bytes
        continue; 
      end
      F2 = fopen( files(n).name , 'rb' );
      if F2 > 0
        F2content = fread( F2 , 'int8=>int8' );
        fclose( F2 );
        if numel(F2content) == numel(Fcontent) && all( F2content == Fcontent )
          dups{ i } = [ dups{i} files(n).name ];
          dups{ n } = [ dups{n} files(i).name ];
          
%           dups{end+1} = files(i).name;
%           break;
        end
      else
        continue;
      end
    end
  end
  fprintf('\n');

end



