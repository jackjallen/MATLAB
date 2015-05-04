function fid = safe_fopen( fname , per,  maxTime , maxTries )

  if nargin < 3 || isempty( maxTime  ), maxTime =  10; end
  
    
  if nargin < 4 || isempty( maxTries ), maxTries = 20; end

  fid = -1;
  
  isf = isfile( fname  );
  
  if strcmp(per,'r') 
    if ~isf
        return;
    end
  end
      
  

  Ntries = 0;
  T0 = clock;
  while etime( clock , T0 ) < maxTime  ||  Ntries < maxTries
    Ntries = Ntries +1;
    
    
    
    fid = fopen( fname , per );
    
    
    if fid > 0,  return; end
    
    T1 = clock; while etime( clock , T1 ) < 0.5, end
  end


end
