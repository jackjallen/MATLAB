function K = FK2Kernel( K , varargin )
  
  K = fftshift( ifftn( K ) );
  
  K = reduceKernel( K , varargin{:} );
  
end
