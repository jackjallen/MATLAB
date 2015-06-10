function CreateSymbolicArray( n , sz )

  n = lower( n );
  N = upper( n );

  if numel(sz) == 1
    inds = (1:sz)';
  else
    szs = {};
    for i = sz
      szs = [ szs 1:i ];
    end
    inds = ndmat( szs{:} );
  end

%   evalin( 'caller' , [ N ' = sym( zeros( [' num2str(sz) '] ));'] )

  for i = 1:size( inds , 1 )
    ind = inds(i,:);
    evalin('caller' , [ N '(' ...
      sprintf( '%d,', ind(1:end) ) '1' ...
      ') = sym( ''' ...
      n  ...
      sprintf('%d',ind) ...
      ''' );' ] );
  end

end
