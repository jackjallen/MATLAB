function x = DCTfilter( x ,  f )

  try
    x = frtn(x , 'dct' );
  catch
    x = dctn(x);
  end

  if iscell( f )
    nd = ndims(x);
    for d = 1:numel( f )
      ff = f{d}; ff = ff(:);

      if size( x , d ) ~= numel(ff)
        error('DCTfilter:invalidFilterSize','the filter have to had the same dim than data.');
      end
      
      ff = ipermute( ff(:) , [ d 1:d-1 d+1:nd ] );
      
      x = bsxfun( @times , x , ff );
      
    end
  else
    x = bsxfun( @times , x , f );
  end

  try    
    x = frtn(x , 'idct' );
  catch
    x = idctn(x);
  end
  
end
