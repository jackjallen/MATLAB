function varargout= argument( a , varargin )
% 
% i= parseargs( a , key1, key2... )       returns the index where are the
%                                         first of the keys
% [a,i]= parseargs( a , key1, key2... )   in i returns the index where are
%                                         the key or 0, and in a the previous a
%                                         without the key
% [a,i,v1..]= parseargs( a , key1, key2... )   
%                                         in i returns the index where are
%                                         any key, and in a the previous a
%                                         without the key, and the follows
%                                         values.
%
% [a,i,v1..]= parseargs( a , key1, key2..., '$DEFS$',v1,v2... )
% 
% [a,i,v1..]= parseargs( a , key1, key2..., {checkfcn,message} )
%

forces = {[],0};
i=1;
while i < numel( varargin )
  try
    if strcmp( varargin{i},'$FORCE$' )
      varargin(i)= [];
      forces= varargin{i};
      varargin(i)= [];
    end
  end
  i= i+1;
end
if numel(forces)==1 || isempty( forces{2} )
  forces{2}=0;
end

for i=1:(nargout-2)
  defaults{i}= [];
end
i=1;
while i < numel( varargin )
  try
    if strcmp( varargin{i},'$DEFS$' )
      for j= 1:numel(defaults)
        varargin(i)= [];
        defaults{j}= varargin{i};
      end
      i= Inf;
    end
  end
  i= i+1;
end


a=a(:);

errormsg = [];
checkfcn = [];

if iscell( varargin{end} )
  checkfcn= varargin{end};
  varargin= varargin(1:end-1);
  if numel( checkfcn ) > 1
    errormsg= checkfcn{2};
  end
  checkfcn= checkfcn{1};
end
v= varargin(:);
zz= 0;

for i=1:numel(a)
  for j=1:numel(v)
    if ischar(a{i}) && ischar(v{j}) && strcmpi( a{i}, v{j} )
      if ~zz
        zz= i;
      end
    end
  end
end


if nargout <= 1 
  if zz > 0
    if isempty( forces{1} )
      varargout{1}= zz;
    else
      varargout{1}= forces{1};
    end
  else
    varargout{1}= forces{2};
  end
else
  varargout{1}= a';

  if zz > 0
    if isempty( forces{1} )
      varargout{2}= zz;
    else
      varargout{2}= forces{1};
    end
  else
    varargout{2}= forces{2};
  end
  
  for i=1:nargout-2
    varargout{i+2}=defaults{i};
  end
  if zz > 0
    a(zz)= [];
    for i=1:nargout-2
      if ~isempty( checkfcn )
        if ~checkfcn( a{zz} )
          if isempty( errormsg )
            error( 'Bad argument parsing at %d-th arg. Check the sintaxis.', zz+i );
          else
            error( ['Bad argument parsing at %d-th arg. ' errormsg ] , zz+i );
          end
        end
      end
      varargout{i+2} = a{zz};
      a(zz)= [];
    end
    varargout{1}= a';
  end
end

