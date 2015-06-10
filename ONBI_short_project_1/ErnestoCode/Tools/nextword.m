function w = nextword( w , letters )

  if nargin < 2
    letters = 'abcdefghijklmnopqrstuvwxyz0123456789_';
  end

  if numel( w ) == 0
    w = letters(1);
  elseif w(end) == letters(end)
    w = [ nextword(w(1:end-1),letters)  letters(1) ];
  else
    w(end) = letters( find(letters==w(end),1) + 1);
  end
  
  
%   if w(end) < 'z' && w(end) >= 'a'
%     w(end) = w(end)+1;
%   elseif w(end) == 'z'
%     w(end) = '0';
%   elseif w(end) >= '0' && w(end) < '9'
%     w(end) = w(end)+1;
%   elseif w(end) == '9'
%     w(end) = '_';
%   elseif w(end) == '_' && numel(w) > 1
%     w = [ nextword(w(1:end-1)) 'a' ];
%   else
%     w = 'aa';
%   end



if 0
%%  
%   clc
%   fclose( fopen( 'feature.txt' , 'w' ) );
  w = 'hbksk_'; 
  while 1
    w=nextword(w); 
    try
      if rand(1)<0.00001
        fprintf('            :::::%s\n',w);
        drawnow;
        fid = fopen('feature.txt','a'); fprintf(fid,'            :::::%s\n',w); fclose(fid);
      end

      feature(w);
      fprintf('%s\n',w);
      fid = fopen('feature.txt','a'); fprintf(fid,'%s\n',w); fclose(fid);
    end
  end  
%%
end


if 0
%%
%   fclose( fopen('commandos.txt','w') );
%   clc
  w = 'ca99p';
  while 1
    w = nextword(w);
    try

      if rand(1)<0.0001
        fprintf('            :::::%s\n',w);
        drawnow;
        fid = fopen('commandos.txt','a');
        fprintf(fid,'            :::::%s\n',w);
        fclose(fid);
      end
      
      if w(1) >= '0' && w(1) <= '9' , continue; end
      if all( w >= '0' ) && all( w <= '9' ), continue; end
      if w(1) == '_', continue; end
      d = which(w);
      if ~isempty(regexp(d,'undoc'))
        fprintf('%s   -> %s\n',w,d);
        fid = fopen('commandos.txt','a');
        fprintf(fid,'%s   -> %s\n',w,d);
        fclose(fid);
        continue;
      end      
      if isempty(d), continue; end
      if strcmpi(d(end-1:end),'.m'), continue; end
      if strncmpi(d,'g:\',3), continue; end
      if strncmpi(d,'/extra/',7), continue; end

      if strcmp(w,'exit'), continue; end
      if strcmp(w,'quit'), continue; end
      if strcmp(w,'clc') , continue; end
      if strcmp(w,'pause') , continue; end
      if strcmp(w,'doc') , continue; end
      if strcmp(w,'hist'), continue; end
      if strcmp(w,'clear'), continue; end

      try
        evalc(w);
        fprintf('%s   -> %s\n',w,d);
        fid = fopen('commandos.txt','a');
        fprintf(fid,'%s   -> %s\n',w,d);
        fclose(fid);
      catch LE
        switch LE.identifier
          case 'MATLAB:UndefinedFunction',
          case 'MATLAB:m_illegal_character',
          otherwise
            fid = fopen('commandos.txt','a');
            fprintf(fid,'%s   -> %s                       %s   --- %s\n',w,d,LE.identifier,LE.message );
            fclose(fid);
        end
      end
    catch
      fid = fopen('commandos.txt','a');
      fprintf(fid,'- %s\n',w );
      fclose(fid);
    end
    
  end  
%%
end



end