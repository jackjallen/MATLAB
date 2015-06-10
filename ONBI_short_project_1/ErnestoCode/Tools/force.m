function c= force( a , op , b ,varargin )

  switch op
    case '+' , op= @plus;
    case '-' , op= @minus;
    case '.*', op= @times;
    case '*' , op= @mtimes;
    case './', op= @rdivide;
    case '.\', op= @ldivide;
    case '/' , op= @mrdivide;
    case '\' , op= @mldivide;
    case '.^', op= @power;
    case '^' , op= @mpower;
    case '<' , op= @lt;
    case '>' , op= @gt;
    case '<=', op= @le;
    case '>=', op= @ge;
    case '~=', op= @ne;
    case '==', op= @eq;
    case '&' , op= @and;
    case '|' , op= @or;
  end

  if ischar(a)
    aa=a;
    try, a= str2num(a); end
    if isempty(a), a= evalin('caller',aa); end
  end
  if ischar(b)
    bb=b;
    try, b= str2num(b); end
    if isempty(b), b= evalin('caller',bb); end
  end

  dims = max( ndims(a) , ndims(b) );
  for d=1:dims, sza(d)= size(a,d); szb(d)= size(b,d); end

  [varargin,i,szc] = parseargs(varargin,'Size','sz','$DEFS$','first');
  switch szc
    case 'first' , szc= sza;
    case 'second', szc= szb;
    case 'max'   , szc= max(sza,szb);
    case 'min'   , szc= min(sza,szb);
  end
  
  [varargin,i,pad] = parseargs(varargin,'Padding','pad','$DEFS$','repmat');
  
  for d=1:dims
    if sza(d) < szc(d)
      a= padding( a , d , szc(d)-sza(d) , pad );
    elseif sza(d) > szc(d)
      a= permute(a,[d setdiff(1:dims,d)]);
      a= a( 1:szc(d) ,:,:,:,:,:,:,:,:,:,:,:,:,:,:,:,:,:,:,:,:,:,:,:,:,:,:,:);
      a= ipermute(a,[d setdiff(1:dims,d)]);
    end
    if szb(d) < szc(d)
      b= padding( b , d , szc(d)-szb(d) , pad );
    elseif szb(d) > szc(d)
      b= permute(b,[d setdiff(1:dims,d)]);
      b= b( 1:szc(d) ,:,:,:,:,:,:,:,:,:,:,:,:,:,:,:,:,:,:,:,:,:,:,:,:,:,:,:);
      b= ipermute(b,[d setdiff(1:dims,d)]);
    end
  end
  
  c= op(a,b);
end
