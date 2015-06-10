function K= gaussianKernel( varargin )
%
%   K= GaussianKernel( -30:30 , 'sigma',s );
%   K= GaussianKernel( -3:3 , -3:3 , 'sigma',s );
%   K= GaussianKernel( -3:3 , -3:3 , -3:3 , 'sigma',s );
%

  [varargin,i,sigma]= parseargs(varargin,'std','sigma','s','$DEFS$',1);
  [varargin,i,tol  ]= parseargs(varargin,'tol','tolerance','t','$DEFS$',0.05);
  [varargin,normalize]= parseargs(varargin,'norm','normal','normalize','n','$FORCE$', {1,0} );
  [varargin,normalize]= parseargs(varargin,'notnorm','notnormal','notnormalize','notn','nn','$FORCE$', {0,1} );
  
  ndims= numel( varargin );
  if numel(sigma)==1
    sigma= repmat(sigma,1,ndims);
  end
%   sigma(ndims+1)=0; sigma(end)=[];

  varargin = cellfun( @(x) x(:) , varargin , 'UniformOutput',false );
  if ndims == 1
    coords = varargin;
  else
    [ coords{1:ndims} ]= ndgrid( varargin{:} );
  end

  K = coords{1}.^2/sigma(1)^2;
  for d= 2:ndims
    K= K + coords{d}.^2/sigma(d)^2;
  end

  K= exp( - K/2 );
  
  if tol > 0
    tol = max( abs(K(:)) )*tol;
    for d = 1:ndims
      K = permute( K , [ d setdiff( 1:max(ndims,2) , d) ] );
      while all( abs( K( 1 ,:) ) < tol )
        K( 1 ,:) = [];
%         varargin{d}( 1 ) = [];
      end
      while all( abs( K(end,:) ) < tol )
        K(end,:) = [];
%         varargin{d}(end) = [];
      end
      K = ipermute( K , [ d setdiff( 1:max(ndims,2) , d) ] );
    end
  end
  
  if normalize
%     V = mean( diff( varargin{1} ));
%     for d=2:ndims
%       V = V * mean( diff( varargin{d} ));
%     end
    K= K/sum(K(:));
  end

end
