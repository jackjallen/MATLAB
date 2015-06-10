function SPLITS = generaSPLITS( nS , nA , nB )

  n     = nA+nB;
  nPmax = round( exp( sum( log(1:n) ) - sum( log( 1:nA) ) - sum( log(1:nB) ) ) );
  if nA == nB, nPmax = nPmax / 2; end
  
  if isempty( nS ), nS = nPmax; end
  if nPmax < nS
    error('no se pueden hacer mas de %d splits', nPmax );
  end
  
  SPLITS = [ true(1,nA) , false(1,nB) ];
  while size( SPLITS , 1 ) < nS
    [~,nSPLITS] = sort( rand(nS,n) , 2 );
    nSPLITS = nSPLITS <= nA;
    if nA == nB
      r = ~nSPLITS(:,1);
      nSPLITS(r,:) = ~nSPLITS(r,:);
    end
    SPLITS = [ SPLITS ; nSPLITS ];
  
    [~,r] = unique( SPLITS , 'rows' , 'first' );
    r = sort(r);
    r = r(1:min(end,nS));
    SPLITS = SPLITS(r,:);
  end
  
  SPLITS = SPLITS.';

end

% 
% 
% 
%   %%    
%   %%    
%   %%    
%   %%    
%   %%    
%   %%    
%   %%    
%   %%    
%   %%    
%   %%    
%   %%    
%   %%    
%   %%    
%   %%    
%   %%    
%   %%    
%   %%    
%   %%    
%   %%    
%   %%    
% 
%   
% nA = 10; nB=5; nS = 250;
% n = nA + nB;
% G  = [ true(1,nA) , false(1,nB) ];
% [~,SPLITS]=sort(rand(nS,n),2);
% SPLITS = G(SPLITS);
% 
% tic
% B = cat(2, SPLITS , false( nS , 32*ceil(n/32) - n ) );
% B = B.';
% B = bin2num( B );
% B = typecast( B , 'int32' );
% B = reshape( B , [] , nS );
% B = B.';
% [~,B]=unique( B , 'rows' ,'stable');
% uniqueSPLITS = SPLITS( B , : );
% toc
% 
% tic
% uniqueSPLITS2 = unique(SPLITS,'rows','stable');
% toc
% 
% isidentical( uniqueSPLITS , uniqueSPLITS2 )
% 
%     %%      
