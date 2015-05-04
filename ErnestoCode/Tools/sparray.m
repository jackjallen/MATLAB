classdef sparray
  %sparray - A class of N-dimensional sparse arrays.
  %
  % USAGE:
  %
  %   S = sparray(X) where X is an ordinary MATLAB sparse matrix converts X into
  %     an "sparray" object. S can be reshaped into an N-dimensional sparse array using
  %     its RESHAPE method, for arbitrary N.
  %
  %   S = sparray( X , [M,N,P,...] ) is equivalent to reshape( sparray(X) , [M,N,P,...] ).
  %
  %

  properties (Access=private)
    data;
    sz;
  end

  methods ( Hidden = true )
    %sparray class constructor
    function obj = sparray( data , varargin )
      %sparray class constructor
      %
      % OBJ=sparray(A) where A is an ordinary MATLAB sparse matrix converts A into
      % an sparray object. OBJ can be reshaped into an N-dimensional sparse array
      % using its RESHAPE method, for arbitrary N.
      %
      % OBJ=sparray(A,[M,N,P,...]) is equivalent to reshape(sparray(A),[M,N,P,...])
      % OBJ=sparray(A, M,N,P,... ) is equivalent to reshape(sparray(A), M,N,P,... )
      % OBJ=sparray(A, M,[],P,...) is equivalent to reshape(sparray(A), M,[],P,...])


      if nargin == 0
        data = 0;
        siz = 1;
      elseif nargin == 1
        siz = size( data );
      elseif nargin == 2
        siz = varargin{1};
      else
        emp = find( cellfun('isempty',varargin) );
        switch numel( emp )
          case 0
            siz = [ varargin{:} ];
          case 1
            siz = [ varargin{ 1:(emp-1) } , 1 , varargin{ (emp+1):end } ];
            missing = numel( data )/prod( siz );
            if mod(missing,1), error 'Bad reshape args'; end
            siz(emp) = missing;
          otherwise
            error 'Too many empty args'
        end
      end

      if numel( siz ) == 0, error('incorrect constructor'); end
      
      siz = fix_size( siz );
      obj.sz   = siz;

      if numel( data ) == 0
        siz = size2D( siz );
        obj.data = sparse([],[],[],siz(1),siz(2));
        return;
      end

      if prod( siz ) ~= numel(data)
        error('The number of DATA elements is not consistent with the given sz.');
      end

      if isa( data , 'sparray' )
        obj.data = reshape( data.data , size2D( siz ) );
        return
      end

      if isequal( size(data) , size2D(siz) )
        obj.data = data;
      else
        obj.data = reshape( data , size2D( siz ) );
      end
      
      if ~issparse(data),
        if ~isa( obj.data ,'double' ) && ~isa( obj.data , 'logical' )
          obj.data = double( obj.data );
        end

        siz = size( obj.data );
        IS = find( obj.data );
        IS = ind2subv( siz , IS );
        try
          VS = obj.data( ~~obj.data );
        catch
          VS = builtin('nonzeros', obj.data );
        end
        
        obj.data = sparse( IS(:,1) , IS(:,2) , VS , siz(1) , siz(2) );
      end

    end

    function varargout = size( obj , dim ) %to hidden
      siz = obj.sz;

      if nargout <= 1 && nargin > 1
        lastDIM = max(dim(:));
        siz( end+1:lastDIM ) = 1;
        varargout{1} = siz( dim );
      elseif nargout <= 1
        varargout{1} = siz;
      elseif nargin > 1
        error('incorrect calling to size');
      else
        siz( end+1:nargout ) = 1;
        siz = [ siz(1:nargout-1) , prod( siz(nargout:end) ) ];
        varargout = num2cell( siz );
      end
    end

    function obj = reshape( obj , varargin )
      
      obj = reshape0( obj , varargin{:} );
      
      siz = size2D( obj.sz );
      
      if ~isequal( size(obj.data) , siz )
        obj.data = reshape( obj.data , siz );
      end
    end
    function obj = reshape0( obj , varargin )
      if nargin == 1, error('incorrect!!'); end
      if nargin == 2
        siz = varargin{1};
      else
        emp = find( cellfun('isempty',varargin) );
        switch numel( emp )
          case 0
            siz = [ varargin{:} ];
          case 1
            siz = [ varargin{ 1:(emp-1) } , 1 , varargin{ (emp+1):end } ];
            missing = numel( obj.data )/prod( siz );
            if mod(missing,1), error 'Bad reshape args'; end
            siz(emp) = missing;
          otherwise
            error 'Too many empty args'
        end
      end
      siz = fix_size( siz );
      if prod( siz ) ~= numel( obj.data )
        error('incorrect sizes' );
      end
      
      obj.sz   = siz;
    end
    
    function obj = vec( obj , d )
      if nargin < 2, d = 1; end
      obj.sz = fix_size( [ ones(1,d-1) , numel(obj.data) ] );
      obj.data = reshape( obj.data , size2D( obj.sz ) );
    end

    function out = full( obj )
      try
        out = reshape( full(obj.data) , obj.sz );
      catch
        error('too large???');
      end
    end

    function out = isidentical( L , R )
      out = false;
      if ~isa(L,'sparray'),   return; end
      if ~isa(R,'sparray'),   return; end
      if ~isidentical( L.sz , R.sz ), return; end
      if ~isidentical( L.data , R.data ), return; end
      out = true;
    end
    
    
    function out = isvector( obj )
      out = numel( obj.sz ) == 2  &&  ( obj.sz(1) == 1 || obj.sz(2) == 1 );
    end
    function out = isscalar(obj)    , out = isscalar(obj.data);         end
    function out = islogical(obj)   , out = islogical(obj.data);        end
    function out = issparse(obj)    , out = issparse(obj.data);         end
    function out = isnumeric(obj)   , out = isnumeric(obj.data);        end
    function out = isfloat(obj)     , out = isfloat(obj.data);          end
    function out = isreal(obj)      , out = isreal(obj.data);           end
    function out = numel(obj)       , out = numel(obj.data);            end
    function out = length( obj )    , out = max( obj.sz );              end
    function out = ndims(obj)       , out = numel( obj.sz );            end
    function out = nnz(obj)         , out = nnz(obj.data);              end
    function out = nzmax(obj)       , out = nzmax(obj.data);            end
    function out = nonzeros(obj)    , out = nonzeros(obj.data);         end
    function obj = real(obj)        , obj.data = real(obj.data);        end
    function obj = imag(obj)        , obj.data = imag(obj.data);        end
    function obj = conj(obj)        , obj.data = conj(obj.data);        end
    function obj = abs(obj)         , obj.data = abs(obj.data);         end
    function obj = sqrt(obj)        , obj.data = sqrt(obj.data);        end
    function obj = logical(obj)     , obj.data = logical(obj.data);     end
    function obj = isnan(obj)       , obj.data = isnan(obj.data);       end
    function obj = isinf(obj)       , obj.data = isinf(obj.data);       end
    function obj = isfinite(obj)    , obj.data = isfinite(obj.data);    end
    function obj = spones(obj)      , obj.data = spones(obj.data);      end
    function obj = double(obj)      , obj.data = double(obj.data);      end
    function obj = uplus(obj)       , obj.data = uplus( obj.data );     end
    function obj = uminus(obj)      , obj.data = uminus( obj.data );    end
    function obj = not(obj)         , obj.data = not( obj.data );       end

    function obj = sign(obj)        , obj.data = sign( obj.data );      end
    function obj = exp(obj)         , obj.data = exp( obj.data );       end
    function obj = expm1(obj)       , obj.data = expm1( obj.data );     end
    function obj = pow2(obj)        , obj.data = pow2( obj.data );      end
    function obj = log(obj)         , obj.data = log( obj.data );       end
    function obj = log2(obj)        , obj.data = log2( obj.data );      end
    function obj = log10(obj)       , obj.data = log10( obj.data );     end
    function obj = reallog(obj)     , obj.data = reallog( obj.data );   end
    function obj = realsqrt(obj)    , obj.data = realsqrt( obj.data );  end
    function obj = cos(obj)         , obj.data = cos( obj.data );       end
    function obj = sin(obj)         , obj.data = sin( obj.data );       end
    function obj = tan(obj)         , obj.data = tan( obj.data );       end
    function obj = cosh(obj)        , obj.data = cosh( obj.data );      end
    function obj = sinh(obj)        , obj.data = sinh( obj.data );      end
    function obj = tanh(obj)        , obj.data = tanh( obj.data );      end
    function obj = acos(obj)        , obj.data = acos( obj.data );      end
    function obj = asin(obj)        , obj.data = asin( obj.data );      end
    function obj = atan(obj)        , obj.data = atan( obj.data );      end
    function obj = acosh(obj)       , obj.data = acosh( obj.data );     end
    function obj = asinh(obj)       , obj.data = asinh( obj.data );     end
    function obj = atanh(obj)       , obj.data = atanh( obj.data );     end
    function obj = round(obj)       , obj.data = round( obj.data );     end
    function obj = fix(obj)         , obj.data = fix( obj.data );       end
    function obj = ceil(obj)        , obj.data = ceil( obj.data );      end
    function obj = floor(obj)       , obj.data = floor( obj.data );     end
    
    
    function obj = mod(L,R)         , obj = bin_op(L,R, @mod );         end
    function obj = rem(L,R)         , obj = bin_op(L,R, @rem );         end
    function obj = plus(L,R)        , obj = bin_op(L,R, @plus );        end
    function obj = minus(L,R)       , obj = bin_op(L,R, @minus );       end
    
%     function L = times(L,R)       , 
%       if ~isa(R,'sparray')  &&  numel(R) == 1
%         L.data = times( L.data , R );
%       else, L = bin_op(L,R, @times );
%       end
%     end
    function obj = times(L,R)       , obj = bin_op(L,R, @times );       end

%     function L = rdivide(L,R)       , 
%       if ~isa(R,'sparray')  &&  numel(R) == 1
%         L.data = rdivide( L.data , R );
%       else, L = bin_op(L,R, @rdivide );
%       end
%     end
    function obj = rdivide(L,R)     , obj = bin_op(L,R, @rdivide );     end
    
    function obj = ldivide(L,R)     , obj = bin_op(L,R, @ldivide );     end
    function obj = power(L,R)       , obj = bin_op(L,R, @power );       end
    function obj = realpow(L,R)     , obj = bin_op(L,R, @realpow );     end
    function obj = eq(L,R)          , obj = bin_op(L,R, @eq );          end
    function obj = ne(L,R)          , obj = bin_op(L,R, @ne );          end
    function obj = gt(L,R)          , obj = bin_op(L,R, @gt );          end
    function obj = ge(L,R)          , obj = bin_op(L,R, @ge );          end
    function obj = lt(L,R)          , obj = bin_op(L,R, @lt );          end
    function obj = le(L,R)          , obj = bin_op(L,R, @le );          end
    function obj = and(L,R)         , obj = bin_op(L,R, @and );         end
    function obj = or(L,R)          , obj = bin_op(L,R, @or );          end
    function obj = xor(L,R)         , obj = bin_op(L,R, @xor );         end


    function out = sparse(obj)      , out = obj.data; end
    function obj = spfun(fun,obj)   , obj.data = spfun(fun,obj.data); end

    function varargout = unique( obj , varargin )
      [varargout{1:nargout}] = unique( obj.data , varargin{:} );
    end

    %%Methods for 2D Matrices
    function spy( obj )
      if numel(obj.sz)~=2, error('Operation spy defined for 2D arrays only'); end
      spy( obj.data );
    end
    function obj = det( obj )
      if numel(obj.sz)~=2, error('Operation det defined for 2D arrays only'); end
      obj = det( obj.data );
    end
    function obj = trace( obj )
      if numel(obj.sz)~=2, error('Operation trace defined for 2D arrays only'); end
      obj = full( trace( obj.data ) );
    end
    function obj = inv( obj )
      if numel(obj.sz)~=2, error('Operation inv defined for 2D arrays only'); end
      obj.data = inv( obj.data );
    end
    function obj = pinv( obj )
      if numel(obj.sz)~=2, error('Operation pinv defined for 2D arrays only'); end
      obj.data = pinv( obj.data );
    end
    function obj = triu( obj , varargin )
      if numel(obj.sz)~=2, error('Operation triu defined for 2D arrays only'); end
      obj.data = triu( obj.data , varargin{:} );
    end
    function obj = tril( obj , varargin )
      if numel(obj.sz)~=2, error('Operation tril defined for 2D arrays only'); end
      obj.data = tril( obj.data , varargin{:} );
    end
    function obj = transpose( obj )
      if numel(obj.sz)~=2, error('Operation transpose defined for 2D arrays only'); end
      obj.data = transpose( obj.data );
    end
    function obj = ctranspose( obj )
      if numel(obj.sz)~=2, error('Operation ctranspose defined for 2D arrays only'); end
      obj.data = ctranspose( obj.data );
    end
    function obj = rot90( obj )
      if numel(obj.sz)~=2, error('Operation rot90 defined for 2D arrays only'); end
      obj.data = rot90( obj.data );
    end

    
    function obj = kron( L , R )
      if       isa( L , 'sparray')  &&  ~isa( R , 'sparray')

        if numel( L.sz ) ~= 2 ||  ndims( R ) ~= 2
          error('Operation kron defined for 2D arrays only');
        else
          obj = sparray( kron( L.data , R ) );
        end

      elseif  ~isa( L , 'sparray')  &&   isa( R , 'sparray')

        if ndims( L ) ~= 2 ||  numel( R.sz ) ~= 2
          error('Operation kron defined for 2D arrays only');
        else
          obj = sparray( kron( L , R.data ) );
        end

      elseif   isa( L , 'sparray')  &&  isa( R , 'sparray')

        if numel( L.sz ) ~= 2 ||  numel( R.sz ) ~= 2
          error('Operation kron defined for 2D arrays only');
        else
          obj = sparray( kron( L.data , R.data ) );
        end

      end
    end

    function obj = mtimes( L , R )
      if       isa( L , 'sparray')  &&  ~isa( R , 'sparray')

        if      numel(R) == 1
          %sparray(rand(3,4))*5
          obj = L;
          obj.data = obj.data * R;
        elseif  numel(L.data) == 1
          %sparray(5)*rand(3,4)
          obj = sparray( L.data * R );
        elseif numel( L.sz ) ~= 2 ||  ndims( R ) ~= 2
          error('Operation mtimes defined for 2D arrays only');
        elseif L.sz(2) == size(R,1)
          %sparray(rand(3,4))*rand(4,5)
          try
            obj = L.data * R;
          catch
            obj = sparse_mtimes( L.data , R , L.sz , size(R) , size2D([ L.sz(1) , size(R,2) ]) );
          end
          if isnumeric( obj ), obj = sparray( obj , [ L.sz(1) , size(R,2) ] ); end
        else
          %sparray(rand(3,4))*rand(5,4)
          error('Inner dimensions must agree.');
        end

      elseif  ~isa( L , 'sparray')  &&   isa( R , 'sparray')

        if      numel(R.data) == 1
          %rand(3,4)*sparray(5)
          obj = sparray( L * R.data );
        elseif  numel(L) == 1
          %5*sparray(rand(3,4))
          obj = R;
          obj.data = L * obj.data;
        elseif ndims( L ) ~= 2 ||  numel( R.sz ) ~= 2
          error('Operation mtimes defined for 2D arrays only');
        elseif size(L,2)  == R.sz(1)
          %sparray(rand(3,4))*rand(4,5)
          obj = sparray( L * R.data , [ size(L,1) , R.sz(2) ] );
        else
          %sparray(rand(3,4))*rand(5,4)
          error('Inner dimensions must agree.');
        end

      elseif   isa( L , 'sparray')  &&  isa( R , 'sparray')

        if      numel(R.data) == 1
          %sparray(rand(3,4))*sparray(5)
          obj = L;
          obj.data = obj.data * R.data;
        elseif  numel(L.data) == 1
          %sparray(5)*sparray(rand(3,4))
          obj = R;
          obj.data = L.data * obj.data;
        elseif numel( L.sz ) ~= 2 ||  numel( R.sz ) ~= 2
          error('Operation mtimes defined for 2D arrays only');
        elseif L.sz(2) == R.sz(1)
          %sparray(rand(3,4))*sparray(rand(4,5))
%           try
%             obj = sparray( L.data * R.data , [ L.sz(1) , R.sz(2) ] );
%           catch
            obj = sparray( sparse_mtimes( L.data , R.data , L.sz , R.sz , size2D( [ L.sz(1) , R.sz(2) ] ) ) , [ L.sz(1) , R.sz(2) ] );
%           end
        else
          %sparray(rand(3,4))*sparray(rand(5,4))
          error('Inner dimensions must agree.');
        end

      end
    end
    
    function obj = mrdivide( L , R )
      if       isa( L , 'sparray')  &&  ~isa( R , 'sparray')

        if      numel(R) == 1
          %sparray(rand(3,4))/5
          obj = L;
          obj.data = mrdivide( obj.data , R );
        else
          error('dividir solamente por escalares!!!');
        end

      elseif  ~isa( L , 'sparray')  &&   isa( R , 'sparray')

        if      numel(R.data) == 1
          %rand(3,4)/sparray(5)
          obj = sparray( mrdivide( L , R.data ) );
        else
          error('dividir solamente por escalares!!!');
        end

      elseif   isa( L , 'sparray')  &&  isa( R , 'sparray')

        if      numel(R.data) == 1
          %sparray(rand(3,4))/sparray(5)
          obj = L;
          obj.data = mrdivide( obj.data , R.data );
        else
          error('dividir solamente por escalares!!!');
        end

      end
    end

    function obj = mpower( L , R )
      if       isa( L , 'sparray')  &&  ~isa( R , 'sparray')
        if numel(L.sz) ~= 2, error('Operation mpower defined for 2D arrays only'); end
        if L.sz(1) ~= L.sz(2), error('sparray must be square'); end
        if  numel(R) == 1
          %sparray(rand(3,3))^2
          obj = L;
          obj.data = mpower( obj.data , R );
        else
          error('mpower solamente para exponentes escalares!!!');
        end

      elseif  ~isa( L , 'sparray')  &&  isa( R , 'sparray')
        if ndims(L) ~= 2, error('Operation mpower defined for 2D arrays only'); end
        if size(L,1) ~= size(L,2), error('sparray must be square'); end

        if      numel(R.data) == 1
          obj = sparray( mpower( L , R.data ) );
        else
          error('mpower solamente para exponentes escalares!!!');
        end

      elseif   isa( L , 'sparray')  &&  isa( R , 'sparray')
        if numel(L.sz) ~= 2, error('Operation mpower defined for 2D arrays only'); end
        if L.sz(1) ~= L.sz(2), error('sparray must be square'); end

        if      numel(R.data) == 1
          obj = L;
          obj.data = mpower( obj.data , R.data );
        else
          error('mpower solamente para exponentes escalares!!!');
        end

      end
    end

    function obj = permute( obj , p )

      siz = obj.sz;
      if ~isequal( sort(p) , 1:max(p) ), error('bad permutation order'); end
      if max(p) < numel( obj.sz ),       error('bad permutation order'); end

      if isequal( p , 1:numel(p) ), return; end
      
      siz(end+1:numel(p)) = 1;

      [C,VS] = get_Coordinates_Values( obj , siz );

      C = C(:,p);
      siz = siz(p);
      siz = fix_size( siz );

      obj = build_SPARRAY( VS , C , siz );
    end

    function obj = circshift( obj , c )
      siz = obj.sz;
      c(end+1:numel(siz)) = 0;
      c = c( 1:numel(siz));
      
      [C,VS] = get_Coordinates_Values( obj , siz );
      for d = 1:numel(c)
        if c(d) == 0, continue; end
        
        cinds = circshift( 1:obj.sz(d) , [ 0 , -c(d) ] );
        C(:,d) = cinds( C(:,d) );
      end
      
      obj = build_SPARRAY( VS , C , siz );
    end

    function obj = ipermute( obj , iperm )
      [null,perm] = sort(iperm);
      obj = permute( obj , perm );
    end

    function obj = squeeze( obj )
      if numel( obj.sz ) <= 2, return; end
      siz = obj.sz;
      siz( siz == 1 ) = [];
      obj.sz = fix_size(siz);
    end

    function obj = flipdim( obj , d )

      if numel(d) ~= 1, error('d expected to be a dimension'); end
      if d<1 || mod(d,1), error('d expected to be a dimension'); end

      siz = obj.sz;
      
      if d > numel( siz )  ||  siz(d) == 1
        return;
      end

      [C,VS] = get_Coordinates_Values( obj );
      
      C(:,d) = obj.sz(d) - C(:,d) + 1;

      obj = build_SPARRAY( VS , C , siz );
    end
    function obj = flipdims( obj , dims )
      if nargin < 2, dims = 1:numel( obj.sz ); end

      if any(dims<1) || any( mod(dims,1) ), error('dims expected to be dimensions'); end

      siz = obj.sz;
      dims = dims( dims <= numel(siz) );

      if isempty( dims ), return; end;
      
      [C,VS] = get_Coordinates_Values( obj );
      
      for d = dims
        if siz(d) > 1
          C(:,d) = siz(d) - C(:,d) + 1;
        end
      end

      obj = build_SPARRAY( VS , C , siz );
    end

    function obj = flipud(obj)
      if numel( obj.sz )~= 2, error('X must be a 2-D matrix.'); end
      obj = flipdim( obj , 1 );
    end
    function obj = fliplr(obj)
      if numel( obj.sz )~= 2, error('X must be a 2-D matrix.'); end
      obj = flipdim( obj , 2 );
    end

    function SS = uneval( obj )
      if nnz(obj.data)*3 > numel(obj)
        S_data = uneval( full(obj.data) );
      else
        S_data = uneval( obj.data );
      end
      
      if isequal( size(obj.data) , obj.sz )
        S = sprintf('sparray(%s)',S_data);
      else
        S_sz   = uneval( obj.sz   );
        S = sprintf('sparray(%s,%s)',S_data,S_sz);
      end

      if nargout > 0
        SS = S;
      else
        Iname = inputname(1);
        if isempty( Iname  ),  fprintf('%s\n', S );
        else,                  fprintf('%s = %s\n', Iname , S );
        end
      end
    end

    function out = isequal( A , varargin )
      out = false;
      
      for i = 1:numel(varargin)
        B = varargin{i};
        
        if        isa( A , 'sparray' )  &&  isa( B , 'sparray')
          thiseq = isequal( A.sz , B.sz ) &&...
                   isequal( reshape( A.data , size2D( A.sz ) ) ,...
                            reshape( B.data , size2D( B.sz ) ) );
        elseif   ~isa( A , 'sparray' )  &&  isa( B , 'sparray')
          thiseq = isequal( size(A) , B.sz ) &&...
                   isequal( reshape( A      , size2D( size(A) ) ) ,...
                            reshape( B.data , size2D( B.sz    ) ) );
        elseif    isa( A , 'sparray' )  && ~isa( B , 'sparray')
          thiseq = isequal( size(B) , A.sz ) &&...
                   isequal( reshape( B      , size2D( size(B) ) ) ,...
                            reshape( A.data , size2D( A.sz    ) ) );
        elseif   ~isa( A , 'sparray' )  && ~isa( B , 'sparray')
          thiseq = isequal( A , B );
        end
        
        if ~thiseq, return; end
      end
      out = true;
    end

    function out = isequalwithequalnans( A , varargin )
      out = false;
      
      for i = 1:numel(varargin)
        B = varargin{i};
        
        if        isa( A , 'sparray' )  &&  isa( B , 'sparray')
          thiseq = isequal( A.sz , B.sz ) &&...
                   isequalwithequalnans( reshape( A.data , size2D( A.sz ) ) ,...
                            reshape( B.data , size2D( B.sz ) ) );
        elseif   ~isa( A , 'sparray' )  &&  isa( B , 'sparray')
          thiseq = isequal( size(A) , B.sz ) &&...
                   isequalwithequalnans( reshape( A      , size2D( size(A) ) ) ,...
                            reshape( B.data , size2D( B.sz    ) ) );
        elseif    isa( A , 'sparray' )  && ~isa( B , 'sparray')
          thiseq = isequal( size(B) , A.sz ) &&...
                   isequalwithequalnans( reshape( B      , size2D( size(B) ) ) ,...
                            reshape( A.data , size2D( A.sz    ) ) );
        elseif   ~isa( A , 'sparray' )  && ~isa( B , 'sparray')
          thiseq = isequalwithequalnans( A , B );
        end
        
        if ~thiseq, return; end
      end
      out = true;
    end


    function obj = cat( d , obj , varargin )
      if ~isa( obj , 'sparray' )
        obj = sparray( obj );
      end
      siz = obj.sz; siz( end+1:d ) = 1;
      
      [C,VS] = get_Coordinates_Values( obj , siz );
      
      for i = 1:numel( varargin )
        B = varargin{i};
        if ~isa( B , 'sparray' ), B = sparray( B ); end
        szB = B.sz; szB( end+1:d ) = 1;
        if ~isequal( siz( [ 1:(d-1)  (d+1):end ] ) , szB( [ 1:(d-1)  (d+1):end ] ) )
          error('incorrect sizes');
        end
        
        [C_B,VS_B] = get_Coordinates_Values( B , szB );
        
        C_B(:,d) = C_B(:,d) + siz(d);
        C  = [ C  ; C_B  ];
        VS = [ VS ; VS_B ];
        
        siz(d) = siz(d) + szB(d);
      end

      obj = build_SPARRAY( VS , C , siz );
    end
    function obj = horzcat( varargin ), obj = cat(2,varargin{:}); end
    function obj = vertcat( varargin ), obj = cat(1,varargin{:}); end
    function obj = repmat( obj , m , n )
      %se puede hace mejor!!!
      if nargin == 2
        if numel(m)==1
          rep = [m,m];
        else
          rep = m;
        end
      else
        rep = [m,n];
      end
      
      for d = 1:numel(rep)
        A = obj;
        for n = 1:( rep(d) - 1 )
          obj = cat(d,obj,A);
        end
      end
    end

    
    function obj = subsref( obj , S )
      switch S.type
        case '()'
          if      numel( S.subs ) == 0
            obj = obj.data;
          elseif  numel( S.subs ) == 1 && isequal( S.subs{1} , ':' )
            obj = vec( obj );
          elseif  numel( S.subs ) == 2 && isequal( S.subs{1} , ':' ) && isequal( S.subs{2} , ':' )
            obj = reshape( obj , obj.sz(1) , [] );
          elseif  numel( S.subs ) >= numel( obj.sz ) && all(cellfun(@(s) isnumeric(s) && isscalar(s) , S.subs ))
            sub = cellfun( @(x)x , S.subs );
            siz = [ obj.sz , ones(1,numel(sub)-numel(obj.sz)) ];
            if any( sub > siz )
              error('index exceeds object dimensions');
            end
            ind = sub2indv( siz , sub );
            obj = obj.data( ind );
          else
            uneval(S);
            error('subsref por ahora esta reducido de funcionalidad');
          end
        otherwise
          uneval(S);
          error('subsref por ahora esta reducido de funcionalidad, ademas para esta llamada hay que revisar numel');
      end
    end
    
    function obj = subsasgn( obj , S , in )
      if numel( in ) ~= 1
        error('solo funciona para asignar un elemento');
      end

      switch S.type
        case '()'
          if  numel( S.subs ) >= numel( obj.sz ) && all(cellfun(@(s) isnumeric(s) && isscalar(s) , S.subs ))
            sub = cellfun( @(x)x , S.subs );
            siz = [ obj.sz , ones(1,numel(sub)-numel(obj.sz)) ];
            if any( sub > siz )
              [C,VS] = get_Coordinates_Values( obj , siz );
              
              VS = [ VS ; in  ];
              C  = [ C  ; sub ];
              siz = max( siz , sub );
              
              obj = build_SPARRAY( VS , C , siz );
            else
              ind = sub2indv( siz , sub );
              obj.data( ind ) = in;
            end
            
          else
            uneval(S);
            error('bad subsasgn');
          end
          
        otherwise
          uneval(S);
          error('subsasgn no tiene esta funcionalidad, ademas para esta llamada hay que revisar numel');
      end

    end


    
    
    function display(obj)

%       dig = ceil( max( log10( obj.sz ) ) + 1 );
%       dig = sprintf( '%%%dd x ' , dig );

      siz = sprintf( '%d x ', obj.sz);
      siz = siz(1:end-2);

      fprintf('sparray -  size:  [ %s]  (%d dims)\n', siz, numel(obj.sz));

      fprintf('    stored data:  ( %d x %d ) , %s', size( obj.data , 1 ) , size( obj.data , 2 ) , class( obj.data ) );
      if issparse( obj.data )
        fprintf(' (sparse)    nzmax: %d' , nzmax(obj.data) );
      end
      
      fprintf('\n');

      n = nnz( obj.data );
      N = numel( obj.data );
      fprintf('        density:  %g %% (  %d  of  %d )\n', n/N*100 , n , N );
      
      mem = 0;
      switch lower( computer('arch') )
        case {'win32'}
          if islogical( obj.data )
            mem = (12-8+1) * n + 4 * size( obj.data , 2 ) + 4;
          else
            mem = 12 * n + 4 * size( obj.data , 2 ) + 4;
          end
        case {'win64'}
          if islogical( obj.data )
            mem = (16-8+1) * n + 4 * size( obj.data , 2 ) + 4;
          else
            mem = 16 * n + 8 * size( obj.data , 2 ) + 8;
          end
      end
      
      fprintf('    data memory~  %d bytes ( %g Mb )  (%g in full mode)\n', mem , mem/1024/1024 , numel(obj.data)*iff(islogical(obj.data),1,8)/1024/1024 );


      if ~issparse( obj.data )
        fprintf(2,'obj no es sparse!!!\n');
      end
      if ~isequal( size(obj.data) , size2D( obj.sz ) )
        fprintf(2,'revisa los tamanos, porque llegue hasta aca??\n');
      end

      
      if N == 1
        fprintf('\n');
        fprintf('value: ');
        disp( obj.data );
      end

      fprintf('\n');
      
    end
    
    function obj = diag( obj , k )
      if nargin < 2, k = 0; end
      
      if isvector( obj )
        obj.data = diag( obj.data(:) , k );
        obj.sz   = size( obj.data );
      else
        if numel( obj.sz ) ~= 2
          error('2D array was expected');
        end
        
        obj.data = diag( obj.data , k );
        obj.sz   = size( obj.data );
        
      end
    end


    function obj = any( obj , d )
      if nargin < 2, d = 1; end
    
      siz = obj.sz;
      siz( end+1:d ) = 1;
      
      rsiz = [ prod(siz(1:d-1)) , siz(d) , prod(siz(d+1:end)) ];
      
      obj.data = logical( obj.data );
      C = get_Coordinates_Values( obj , rsiz );

      C(:,2) = 1;
      siz(d) = 1;
      rsiz(2) = 1;

      C    = ind2subv( siz , unique( sub2indv( rsiz , C ) ) );
      
      obj = build_SPARRAY( true , C , siz );
    end
    
    function obj = all( obj , varargin )     , obj = collapse_op( @all      , obj , varargin{:} ); end
    function obj = sum( obj , varargin )     , obj = collapse_op( @sum      , obj , varargin{:} ); end
    function obj = cumsum( obj , varargin )  , obj = collapse_op( @cumsum   , obj , varargin{:} ); end
    function obj = prod( obj , varargin )    , obj = collapse_op( @prod     , obj , varargin{:} ); end
    function obj = cumprod( obj , varargin ) , obj = collapse_op( @cumprod  , obj , varargin{:} ); end
    function obj = mean( obj , varargin )    , obj = collapse_op( @mean     , obj , varargin{:} ); end
    function obj = median( obj , varargin )  , obj = collapse_op( @median   , obj , varargin{:} ); end

    function obj = diff( obj , n , d )
      if nargin < 2, n = 1; end
      if nargin < 3, d = 1; end
      obj = collapse_op( @diff , obj , d , n );
    end
    function obj = var( obj , w , d )
      if nargin < 2, w = 0; end
      if nargin < 3, d = 1; end
      obj = collapse_op( @var , obj , d , w );
    end
    function obj = std( obj , flag , d )
      if nargin < 2, flag = 0; end
      if nargin < 3, d = 1; end
      obj = collapse_op( @std , obj , d , flag );
    end
    
    function [obj , idx ] = sort( obj , d , varargin )
      if numel( varargin ) == 0        ,     mode = 'ascend';
      elseif numel( varargin ) == 1    ,     mode = lower( varargin{1} );
      else                             ,     error('too many inputs');
      end

      switch mode
        case 'descend'
        case 'ascend'
        otherwise, error('unknown mode');
      end
      
      if nargin < 2, d = 1; end
      if nargout > 1
        [ obj , idx ] = collapse_op( @sort , obj , d );
        dims = 1:max( d , numel(obj.sz ) );
        idx = reshape( idx , obj.sz( dims( [ d , 1:(d-1) , (d+1):end ] ) ) );
        idx = ipermute( idx , dims( [ d , 1:(d-1) , (d+1):end ] ) );
        
        switch mode
          case 'descend',  obj = flipdim( obj , d ); idx = flipdim( idx , d );
          case 'ascend'
          otherwise, error('unknown mode');
        end
        
      else
        obj = collapse_op( @sort , obj , d );

        switch mode
          case 'descend',  obj = flipdim( obj , d );
          case 'ascend'
          otherwise, error('unknown mode');
        end

      end
    end
    
    
    %%metodos que faltan
    function obj = spdiags( obj , k ), error('niy'); end
    function obj = colon( obj ), error('niy'); end
    function obj = bsxfun( obj ), error('niy'); end
    function obj = find( obj ), error('niy'); end

    function obj = eig( obj ), error('niy'); end
    function obj = svd( obj ), error('niy'); end

    function obj = mldivide( obj ), error('niy'); end

    function obj = max( obj ), error('niy'); end
    function obj = min( obj ), error('niy'); end
    function obj = prctile( obj ), error('niy'); end
  
    function obj = trapz( obj ), error('niy'); end
    function obj = cumtrapz( obj ), error('niy'); end

    function obj = blkdiag( obj ), error('niy'); end

    function obj = fftshift( obj ), error('niy'); end
    function obj = ifftshift( obj ), error('niy'); end

    function obj = fftn( obj ), error('niy'); end
    function obj = ifftn( obj ), error('niy'); end
    function obj = fft2( obj ), error('niy'); end
    function obj = ifft2( obj ), error('niy'); end
    function obj = fft( obj ), error('niy'); end
    function obj = ifft( obj ), error('niy'); end
    function obj = dctn( obj ), error('niy'); end
    function obj = idctn( obj ), error('niy'); end
    function obj = dct2( obj ), error('niy'); end
    function obj = idct2( obj ), error('niy'); end
    function obj = dct( obj ), error('niy'); end
    function obj = idct( obj ), error('niy'); end
    function obj = conv( obj ), error('niy'); end
    function obj = conv2( obj ), error('niy'); end
    function obj = convn( obj ), error('niy'); end
    function obj = filter( obj ), error('niy'); end

  
    function obj = intersect( obj ), error('niy'); end
    function obj = ismember( obj ), error('niy'); end
    function obj = issorted( obj ), error('niy'); end
    function obj = setxor( obj ), error('niy'); end
    function obj = union( obj ), error('niy'); end
  
  end
  
  methods ( Access = private , Hidden = true )

    %%Binary ops
    function  L = bin_op(L,R,op)
      if       isa(L,'sparray')   &&  ~isa(R,'sparray')
        if      numel(R) == 1
          %sparray(sprand(10,100,0.1),[10 10 10]) + 10
          L.data = feval(op, L.data , R );
          if ~issparse( L.data ), L.data = sparse(L.data); end
        elseif  numel(L.data) == 1
          %sparray(1) + double( randn([10,20,30])>0.5 )
          L.data = feval(op, L.data , R );
          L.sz   = size( R );
          L.data = reshape( L.data , size2D( L.sz ) );
          L.data = sparse( L.data );
        elseif  ~isequal( L.sz , size(R) )
          error('inconsistent sizes');
        else
          %sparray(1) + 2
          %sparray(sprand(10,600,0.1),[10 20 30]) + double( randn([10,20,30])>0.5 )
          L.data = feval(op, L.data , sparse( reshape( R , size(L.data) ) ) );
        end

      elseif  ~isa(L,'sparray')   &&   isa(R,'sparray')
        if      numel(R.data) == 1
          %double( randn([10,20,30])>0.5 ) + sparray(10)
          R.sz = fix_size( size(L) );
          R.data = feval(op, L , R.data );
          R.data = reshape( R.data , size2D( R.sz ) );
          R.data = sparse( R.data );
          L = R;
        elseif  numel(L) == 1
          %1 + sparray(sprand(10,600,0.1),[10 20 30])
          R.data = feval(op, L , R.data );
          if ~issparse( R.data ), R.data = sparse(R.data); end
          L = R;
        elseif  ~isequal( size(L) , R.sz )
          error('inconsistent sizes');
        else
          %2 + sparray(1)
          %double( randn([10,20,30])>0.5 ) + sparray(sprand(10,600,0.1),[10 20 30])
          R.data = feval(op, sparse( reshape( L , size(R.data) ) ) , R.data );
          L = R;
        end

      elseif   isa(L,'sparray')   &&   isa(R,'sparray')
        if      numel(R.data) == 1
          %sparray(sprand(10,600,0.1),[10 20 30]) + sparray(10)
          L.data = feval(op, L.data , R.data );
        elseif  numel(L.data) == 1
          %sparray(1) + sparray(sprand(10,600,0.1),[10 20 30])
          L.sz = fix_size( R.sz );
          L.data = feval(op, L.data , R.data );
        elseif  ~isequal( L.sz , R.sz )
          error('inconsistent sizes');
        else
          %sparray(2) + sparray(1)
          %sparray(sprand(10,600,0.1),[10 20 30]) + sparray(sprand(10,600,0.1),[10 20 30])
          if isequal( size(L.data) , size(R.data) )
            L.data = feval(op, L.data , R.data );
          else
            warning('algo va mal!!!!');
            L.data = feval(op, reshape( L.data , size2D( L.sz ) ) ,...
              reshape( R.data , size2D( L.sz ) ) );
          end
        end

      end
    end

    
    %%collapse array's dimensions functions
    function [obj, out2] = collapse_op( op , obj , d , varargin )
      if nargin < 2, d = 1; end
      
      dims = 1:max( d , numel(obj.sz) );
      obj = permute( obj , dims( [ d , 1:(d-1) , (d+1):end ] ) );
      siz = obj.sz;
      
      obj.data = reshape( obj.data , siz(1) , prod( siz(2:end) ) );
      
      if nargout > 1
        [ obj.data , out2 ] = feval( op , obj.data , 1 , varargin{:} );
      else
        if isequal( op , @diff )
          obj.data = feval( op , obj.data , varargin{:} , 1 );
        elseif isequal( op , @var )
          obj.data = feval( op , obj.data , varargin{:} , 1 );
        elseif isequal( op , @std )
          obj.data = feval( op , obj.data , varargin{:} , 1 );
        else
          obj.data = feval( op , obj.data , 1 , varargin{:} );
        end
      end

      obj.sz(1) = size(obj.data,1);
      obj.data = reshape( obj.data , size2D( obj.sz ) );
      
      obj = ipermute( obj , dims( [ d , 1:(d-1) , (d+1):end ] ) );
    end
    
    
    function [C,VS] = get_Coordinates_Values( obj , siz )
      if nargin < 2, siz = obj.sz; end

      if nargout > 1
        [I,J,VS] = find( obj.data );
        C = IJ_2_ndCoords( siz , I , J );
      else
        C = ind2subv( siz , find( obj.data ) );
      end
    end
    
  end
 
  methods ( Static = true , Hidden = true )
    
    function out = end( obj , K , N )
      %end - sameindexing rules as for full arrays
      kk = false([ obj.sz , ones( 1 , max(K,N)-numel(obj.sz) ) , 0 ]);
      out = builtin('end', kk , K , N );
    end
    
  end
  
  
end


function obj = build_SPARRAY( VS , C , siz )
  [I,J] = ndCoords_2_IJ( siz , C );
  obj = sparray( sparse( I , J , VS , size2D( siz ,1) , size2D( siz ,2) ) , siz );
end
function C = IJ_2_ndCoords( osz , I , J )
  if nargin == 2
    p = [ 1 ; size2D( osz , 1 ) ];
    I = I - 1;
    I = I*p;
    I = I + 1;
  elseif nargin == 3
    p = size2D( osz , 1 );
    J = J-1;
    J = J*p;
    I = I + J;
  end
  C = ind2subv( osz , I );
end
function [I,J] = ndCoords_2_IJ( osz , C )
  sz = [ osz , ones(1, size(C,2) - numel(osz) )];
  I = sub2indv( sz , C );

  sz = size2D( osz );

%   J  = ceil( I / sz(1) );
%   I  = rem( I , sz(1) );
%   if nargout < 2
%     I = [ I , J ];
%   end

  I = ind2subv( sz , I );
  if nargout > 1
    J = I(:,2);
    I = I(:,1);
  end
end
function siz = fix_size( siz )
  siz = siz( 1:find(siz~=1,1,'last' )     );
  siz = [ siz , ones(1,2-numel(siz)) ];
end
function siz = size2D( siz , d )
  siz = siz(   find(siz~=1,1,'first'):end );
  siz = siz( 1:find(siz~=1,1,'last' )     );
  siz = [ siz , ones(1,2-numel(siz)) ];

%   if nargin == 1
%     siz = [ prod(siz(1:end-1)) , siz(end) ];
%   elseif d == 1
%     siz = prod(siz(1:end-1));
%   elseif d == 2
%     siz = siz(end);
%   else
%     error('????');
%   end


%   siz = cumprod( siz );
%   
%   idx = find( siz(1:end-1) < 2^31 , 1 , 'last' );
%   
%   siz = [ siz(idx) , siz(end)/siz(idx) ];
  

  siz = [ prod(siz(1:end-1)) , siz(end) ];
  
  while siz(1) > 2^31
    f = factor1( siz(1) );
    if f == 1, break; end
    siz = [ siz(1)/f(1) , siz(2)*f(1) ];
  end
  
  if nargin > 1, siz = siz(d); end
  
end









%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function test
if 0

%% sparray
sz = [5 20 10 7 9 3 9 ];

T = sprand(sz(1),prod(sz(2:end)),0.1);
F = reshape( full(T) , sz );

S = sparray( T , sz )
SS = sparray( F , 5 , 20 ,[] , 7 ,9,3,9,1,1,1); isidentical( S , SS )
SS = sparray( F ); isidentical( S , SS )
SS = sparray( F ,sz ); isidentical( S , SS )
SS = reshape( sparray( T ) , sz ); isidentical( S , SS )
SS = reshape( sparray( T ) , sz ); isidentical( S , SS )

N = ndSparse( T ,sz);
isidentical( full(N) , full(S) )

%% size method
clc
for n = 1:15
clear szF szS
[szF{1:n}] = size( F );
[szS{1:n}] = size( S ); uneval(szS)
isidentical( szF , szS )
end

for n = 1:15
szF = size( F , n );
szS = size( S , n )
isidentical( szF , szS )
end

[kk1,kk2] = size( F , 2 )
[kk1,kk2] = size( S , 2 )

%% reshape
reshape( S , 3,4,5,[] )
reshape( S , 3,4,5,[],7,9,10 )
reshape( S , 3,4,[],6,7,9,10 )

isidentical( full(reshape( S , 3,4,5,[],7,9,10 )) , reshape( F , 3,4,5,[],7,9,10 ) )
reshape( S , [] , 1 )
isidentical( full(reshape( S , [],1 )) , F(:) )

isidentical( full( vec( S ) ) , vec( F ) )
isidentical( full( vec( S , 5 ) ) , vec( F , 5 ) )

%%
isscalar( S ),   isscalar( sparray(1) )
islogical( S ),  islogical( ~~S )
issparse( S )   %allways!!!
isnumeric( S ), isnumeric( ~~S )
isfloat( S ), isfloat( ~~S )
ndims( S )
size( nonzeros( S ) ), uneval( nonzeros( ~~S ) )
isidentical( logical(S) , ~~S ), logical(S)
isidentical( full( ~~S ), ~~F )
isidentical( double( ~~S ) , spones( S ) )
isidentical( sparse( S ) , S() )
isidentical( getv(struct(S),'.data'),S() )

%% binary ops, que preservan sizes , por ejemplo +
F = randn([2,3,4,5]); F(F>0) = 0;

% sparray  @op@  numeric scalar
L = sparray( F );     R = 20;
LR = L + R, isidentical( full(LR) , full(L)+full(R) )
% scalar sparray @op@ numeric
L = sparray(20);      R = F;
LR = L + R, isidentical( full(LR) , full(L)+full(R) )
% sparray @op@ numeric
L = sparray( F );     R = F;
LR = L + R, isidentical( full(LR) , full(L)+full(R) )


% numeric  @op@  sparray scalar
L =  F ;     R = sparray(20);
LR = L + R, isidentical( full(LR) , full(L)+full(R) )
% numeric scalar @op@ sparray
L = 20;      R = sparray( F );
LR = L + R, isidentical( full(LR) , full(L)+full(R) )
% numeric @op@ sparray
L = F ;     R = sparray( F );
LR = L + R, isidentical( full(LR) , full(L)+full(R) )


% sparray @op@  sparray scalar
L =  sparray( F );     R = sparray(20);
LR = L + R, isidentical( full(LR) , full(L)+full(R) )
% scalar sparray @op@ sparray
L = sparray( 20 );      R = sparray( F );
LR = L + R, isidentical( full(LR) , full(L)+full(R) )
% sparray @op@ sparray
L = sparray( F );     R = sparray( F );
LR = L + R, isidentical( full(LR) , full(L)+full(R) )


%% Methods for 2D Matrices
F = randn(20,20);
S = sparray( F );
isidentical( det(F) , det(S) )
isidentical( trace(F) , trace(S) )
maxnorm( inv(F) , full( inv(S) ) )
isidentical( triu(F,4) , full( triu(S,4) ) )
isidentical( F' , full( S' ) )

%% mtimes
A = rand(5,4); B = rand(4,3); n = 20;
sA = sparray(A); sB = sparray(B); sn = sparray(n);

% sparray * numeric scalar
sA_n = sA * n, isidentical( full(sA_n) , A*n )
% scalar sparray * numeric
sn_A = sn * A, isidentical( full(sn_A) , n*A )
% sparray * numeric
sA_B = sA * B, isidentical( full(sA_B) , A*B )


% numeric * scalar sparray
A_sn = A * sn, isidentical( full(A_sn) , A*n )
% numeric scalar * sparray
n_sA = n * sA, isidentical( full(n_sA) , n*A )
% numeric * sparray
A_sB = A * sB, isidentical( full(A_sB) , A*B )


% sparray * scalar sparray
sA_sn = sA * sn, isidentical( full(sA_sn) , A*n )
% scalar sparray * sparray
sn_sA = sn * sA, isidentical( full(sn_sA) , n*A )
% sparray * sparray
sA_sB = sA * sB, isidentical( full(sA_sB) , A*B )


%% mrdivide    (numerator: array, denominator: scalar) 
N = randn([10,20,30]); sN = sparray( N );
D = 20; sD = sparray(D);

sN/D , isidentical( full( sN/D  ) , N/D )
 N/sD, isidentical( full(  N/sD ) , N/D )
sN/sD, isidentical( full( sN/sD ) , N/D )

%% permute
F = randn([3,4,1,1,1,2,3,4,2,2,1,2,5,1,1,2]);
S = sparray( F );
perm = shuffle( 1:ndims(S) )

isidentical(  permute(F,perm) , full(  permute(S,perm) ) )
isidentical( ipermute(F,perm) , full( ipermute(S,perm) ) )

for n = 1:ndims(F)+2
isidentical( flipdim(F,n) , full( flipdim(S,n) ) )
end

isequal( circshift(F,[0 2 3 -5 -50]) , circshift(S,[0 2 3 -5 -50]) )


%% squeeze
F = randn([1,1,1,3,4,1,1,1,2,3,4,2,2,1,2,5,1,1,1,1,1,2]); S = sparray(F);
F = randn([1,4]); S = sparray(F);
isidentical( squeeze(F) , full( squeeze(S) ) )


%% cat
F = randn([3,4,5,1,1,6]);
S = sparray(F);

isequal( cat(1,S,S) , cat(1,F,F) )
isequal( cat(1,S,F) , cat(1,F,F) )
isequal( cat(10,S,F) , cat(10,F,F) )

isequal( [[S,S,F];[S,S,S];[F,F,F]] , [[F,F,F];[F,F,F];[F,F,F]] )

isequal( repmat( S , 4 ) , repmat( F , 4 ) )
isequal( repmat( S , 4,3 ) , repmat( F , 4,3 ) )
isequal( repmat( S , 1,4 ) , repmat( F , 1,4 ) )
isequal( repmat( S , [4,1,2,1,5,1,1,1,10,1,1,1,2] ) , repmat( F , [4,1,2,1,5,1,1,1,10,1,1,1,2] ) )

%%
F = randn([3,4,5,1,1,6]);
S = sparray(F);

isequal( S() , getv( struct(S) , '.data' ) )
isequal( S(:) , F(:) )
isequal( S(:,:) , F(:,:) )
isequal( S(1,2,end,1,1,5,1) , F(1,2,end,1,1,5,1) )

S(1,2,3,1,1,5) = pi; F(1,2,3,1,1,5) = pi;
isequal( S , F )

S(1,2,8,1,1,50) = pi; F(1,2,8,1,1,50) = pi;
isequal( S , F )
size( S )


%% "collapsing" functions
F = randn([20,10,8,5]); %F(F<0) = 0;
S = sparray( F );
N = ndSparse( F );

isequal( sum( F , 3 ) , sum( S,3 ) )
isequal( cumsum( F , 3 ) , cumsum( S,3 ) )
isequal( prod( F , 3 ) , prod( S,3 ) )
isequal( cumprod( F , 3 ) , cumprod( S,3 ) )
isequal( mean( F , 3 ) , mean( S,3 ) )
isequal( median( F , 3 ) , median( S,3 ) )


isequal( diff( F ) , diff( S ) )
isequal( diff( F , 3 ) , diff( S , 3 ) )
isequal( diff( F , 4 , 3 ) , diff( S , 4 ,3 ) )


isequal( var( F ) , var( S ) )
isequal( var( F , 0 ) , var( S , 0 ) )
isequal( var( F , 1:size(F,3) , 3 ) , var( S , 1:size(F,3) ,3 ) )


isequal( std( F ) , std( S ) )
isequal( std( F , 0 ) , std( S , 0 ) )
isequal( std( F , 1 , 3 ) , std( S , 1 ,3 ) )



[Ss,Si] = sort( S , 3 , 'descend' );
[Fs,Fi] = sort( F , 3 , 'descend' );

isequal( Ss , Fs ), isequal( Si , Fi )

maxnorm( sum( F , 3 ) , full( sum( N , 3 ) ) )

tic; for kk=1:100,sum(S,3); end; toc
tic; for kk=1:100,sum(N,3); end; toc

tic; all(S,3); toc; tic; all(N,3); toc
tic; any(S,3); toc; tic; any(N,3); toc



%% compare performance !
F = randn([20,30,40,50,10]); F(F<1) = 0;
sz = size(F);
T = sparse( F(:,:) );
S = sparray( F );
N = ndSparse( F );

%%
clc
%       for i=1:5,  FF=F.*100; end;
% tic;  for i=1:20, FF=F.*100; end; toc

      for i=1:5,  TT=T.*100; end;
tic;  for i=1:20, TT=T.*100; end; toc

      for i=1:5,  SS=S.*100; end;
tic;  for i=1:20, SS=S.*100; end; toc

%       for i=1:5,  kk=sparray( T.*100 , sz ); end;
% tic;  for i=1:20, kk=sparray( T.*100 , sz ); end; toc
% 
%       for i=1:5,  NN=N.*100; end;
% tic;  for i=1:20, NN=N.*100; end; toc

isequal( FF , full(SS) )

%%
clc
%       for i=1:5,  FF=F./100; end;
% tic;  for i=1:20, FF=F./100; end; toc

      for i=1:5,  TT=T./100; end;
tic;  for i=1:20, TT=T./100; end; toc

      for i=1:5,  SS=S./100; end;
tic;  for i=1:20, SS=S./100; end; toc

%       for i=1:5,  NN=N./100; end;
% tic;  for i=1:20, NN=N./100; end; toc

isequal( FF , full(SS) )

%%
clc
%       for i=1:5,  FF=F*100; end;
% tic;  for i=1:20, FF=F*100; end; toc

      for i=1:5,  TT=T/100; end;
tic;  for i=1:20, TT=T/100; end; toc

      for i=1:5,  SS=S/100; end;
tic;  for i=1:20, SS=S/100; end; toc

%       for i=1:5,  NN=N*100; end;
% tic;  for i=1:20, NN=N*100; end; toc

isequal( FF , full(SS) )
%%
clc
      for i=1:5,  FF=F+100; end;
tic;  for i=1:20, FF=F+100; end; toc

      for i=1:5,  TT=T+100; end;
tic;  for i=1:20, TT=T+100; end; toc

      for i=1:5,  SS=S+100; end;
tic;  for i=1:20, SS=S+100; end; toc

isequal( FF , full(SS) )


end
end