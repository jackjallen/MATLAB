function x = subarray( x , varargin )

  try

    x = x( varargin{:} );

  catch

    nd = ndims( x );
    inds = repmat( {':'} , [1,nd] );
    for d = 1:numel( varargin )
      idx = varargin{d};

      idx_nans = ( idx < 1 )  |  ( idx > size(x,d) );

      if any( idx_nans )
        idx( idx_nans ) = 1;
        x = x( inds{1:d-1} , idx , inds{d+1:end} );
        x( inds{1:d-1} , idx_nans , inds{d+1:end} ) = NaN;
      else
        x = x( inds{1:d-1} , idx , inds{d+1:end} );
      end
    end

  end

end
