function S = equalfields( a, b , varargin )
  if nargin<1
    S= [];
    return;
  end

  if nargin<2
    S= a;
    return;
  end


  S= [];
  if ~isstruct(a) || ~isstruct(b)
    return;
  end
  
  fnames= fieldnames(a);
  for f= 1:numel(fnames)
    fn = fnames{f};
    if ~isfield(b,fn)
      continue;
    end
    aa= a.(fn);
    bb= b.(fn);
    if ~strcmp( class(aa) , class(bb) )
      continue;
    end
    if isstruct(aa)
      aabb= equalfields(aa,bb);
      if ~isempty( aabb )
        S.(fn)= aabb;
      end
      continue;
    end
    if isequal(aa,bb)
      S.(fn) = aa;
    end
  end

  for i=1:numel(varargin)
    S= equalfields(S,varargin{i});
  end
  
  
end
