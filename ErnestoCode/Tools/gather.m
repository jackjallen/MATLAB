function x = gather(x)

  if isa( x , 'parallel.gpu.GPUArray' )
    x = gather( x );
  end

end
