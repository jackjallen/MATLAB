function [ isf , fname ] = isfile( fname )
% 
% isfile(fname)
% 
% return  1   if fname is an existing file
% return -1   if fname is an existing dir
% return  0   if the path doesn't exist
%   
  
  fname= fixfilename( fname );
  
  d= dir( fname );
  
  if numel(d)==1
      isf= 1;
      return;
  end

  if numel(d) > 1
    isf= -1;
    fname= fixfilename([fname filesep]);
    return;
  end

  if numel(d)==0
    isf= 0;
    while numel( dir(fname) ) <= 1
      fname= fileparts(fname);
      if isempty(fname)
        return;
      end
    end
    fname= fixfilename([fname filesep]);
  end
  
end
