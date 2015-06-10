function L = spectralLaplacian( sz , OPtype )

  if nargin < 2
    OPtype = 'impulse';
  end

  if iscell( sz )
    W  = cellfun( @(x)  2*pi*vec( (1:numel(x)) - ceil(numel(x)/2) )/numel(x)/mean(diff(x(:))) , sz , 'UniformOutput', false );
    sz = cellfun( 'prodofsize' , sz );
  else
    W  = arrayfun( @(n) 2*pi*vec( (1:n) - ceil(n/2) )/n , sz , 'UniformOutput', false );
  end
  W  = ndmat_mx( W{:} );

  L = fftshift( reshape( -sum( W.^2 , 2 ) , [ sz(:).' , 1 ] ) );

  
  switch lower( OPtype )
    case 'impulse'
      L = real( fftshift( ifftn( L ) ) );
    case {'fft','fk'}
      
      
      
  end

end
