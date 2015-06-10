function K = kronsel( A , B , rs , cs )

  if      nargin == 2
    
    K = kron( A , B );

  elseif  nargin == 4

    if islogical( rs ), rs = find( rs ); end
    if islogical( cs ), cs = find( cs ); end
    
    K = kron( A , B );

    if      ~isempty( rs )  &&  ~isempty( cs )
      
%       K = zeros( numel(rs) , numel(cs) );
%       
%       for i = 1:numel(rs)
%         for j = 1:numel(cs)
%           K(i,j) = A( rem(
%         end
%       end
      
      K = K( rs , cs );
    elseif  ~isempty( rs )  &&   isempty( cs )
      K = K( rs , : );
    elseif   isempty( rs )  &&  ~isempty( cs )
      K = K( : , cs );
    elseif   isempty( rs )  &&   isempty( cs )
      K = kron( A , B );
    end
    
    
  else
    
    error('2 or 4 arguments expected');
  end
  
  
end
