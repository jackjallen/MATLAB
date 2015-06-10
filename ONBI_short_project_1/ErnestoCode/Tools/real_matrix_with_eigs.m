function M = real_matrix_with_eigs( L )

  M = {};
  if iscell(L)
    for l = 1:numel(L)
      if isreal( L{l} )
        M{end+1} = L{l};
      else
        M{end+1} = [ real( L{l} ) , imag( L{l} ) ; -imag( L{l} ) , real( L{l} ) ];
      end
    end
    M = blkdiag( M{:} );
  end
  
%   R = expm( skewmatrix( randn( 1 , size(M,1)*(size(M,1)-1)/2 ) ) );
  R = expm( randn(size(M)) );
  
  M = R*M/R;
  


%   if iscell(L)
%     L = cell2mat( L(:) );
%     L = [ L( ~imag(L) ) ; L( ~~imag(L) ) ; vec( L( ~~imag(L) )' ) ];
%   end
%   
%   order = @(e) sortrows( [ real( e(:) ) , imag( e(:) ) ] , [1 2] );
%   L = order( L );
% 
%   edist = @(M) fro2( L - order( eig(M) ) );
%   
%   M = zeros( numel(L)/2 );
%   for reg = [ geospace( 1000 , 1e-10 , 20 ) , 0 ]
%     M = Optimize( @(M) edist( M ) + reg*fro2(M) , M , 'methods',{'conjugate','coordinate',numel(M)},struct('COMPUTE_NUMERICAL_JACOBIAN',{{'5'}}) , struct('COORDINATES_ORDER',@(N) 1:N) );
%   end

end
