function o = seval( expr , varargin )
% 
% Example:
%   A= rand(5);
%   [U,S,V] = svd(A);
%
%   seval svd(A) 
%
%   seval( 'svd(A)' ,'{2}' )
%   seval( 'svd(A)' ,2 )
%
%   S(2,2)
%   seval( 'svd(A)' ,'{2}','(2,2)' )
%   seval( 'svd(A)' ,2,'(2,2)' )
% 
%   V(1,:)
%   seval( 'svd(A)' ,'{3}','(1,:)' )
%   seval( 'svd(A)' ,3,'(1,:)' )
%

  nouts= 0;
  try
    if isnumeric( varargin{1} )
      nouts= varargin{1};
      varargin(1)= [];
    end
  end

  if ~nouts
    nouts= 1;
    while 1
      try
        [o{1:nouts}]= evalin( 'caller', expr );
        nouts= nouts+1;
      catch
        break;
      end
    end
  else
    [o{1:nouts}]= evalin( 'caller', expr );
    o= o{end};
  end

  for i=1:numel(varargin)
    if isnumeric( varargin{i} )
      if iscell( o )
        o = eval( ['o{' num2str( varargin{i} ) '};'] );
      else 
        o = eval( ['o(' num2str( varargin{i} ) ');'] );
      end
    else
      o = eval( ['o' varargin{i}] );
    end
  end
  
%   try
%     for i=1:numel( varargin )
%       if isnumeric( varargin{i} )
%         varargin(i)= { varargin(i) };
%       end
%     end
%     var = subsref( o , substruct( varargin{:} ) );
%   catch
%     var = o;
%   end
end
