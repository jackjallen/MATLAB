function K= sincKernel( varargin )
%
%   K= sincKernel( -30:30 , 'z', z ); (default z = pi/2)
%   K= sincKernel( -3:3 , -3:3 , 'zero', z );
%   K= sincKernel( -3:3 , -3:3 , -3:3 , 'z', z );
%

  [varargin,i,zero]= parseargs(varargin,'Zeros','$DEFS$',pi/2);

  ndims= numel( varargin );
  
  if numel(zero)==1
    zero= repmat(zero,1,ndims);
  end
  zero(end+1)=0; zero(end)=[];

  if ndims == 1
    coords = varargin;
  else
    [ coords{1:ndims} ]= ndgrid( varargin{:} );
  end

  K= sinc( coords{1}/zero(1) );
  for d= 2:ndims
    K= K .* sinc( coords{d}/zero(d) );
  end

end
