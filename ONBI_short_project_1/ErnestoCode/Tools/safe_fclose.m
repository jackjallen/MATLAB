function status_ = safe_fclose( fid )

  status = [];
  
  if ischar( fid )
    
    all_fids = fopen('all');
    
    for f = all_fids
      if isequal( fopen(f) , fid )
        status(end+1) = fclose( f );
      end
    end

  elseif isnumeric( fid )
    
    fid = fopen( fid );
    
    if ~isempty( fid )
      status = safe_fclose( fid );
    end
    
  end

  if nargout
    status_ = status;
  end
  

end

