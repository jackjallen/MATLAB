function [t,p,dfe] = studentsttest2( A , B , varargin )

  dimC = 1;  % each row means a component of the vector
  dimS = 2;  % each subject is a column
  dimT = 3;  % each plane is a different test to be performed

  if ndims(A) ~= ndims(B)        , error('different ndims in A and B.'); end
  if ndims(A) > 3                , error('more than 3 dimensions are not allowed.'); end
  if size(A,dimC) ~= size(B,dimC), error('A and B must have the same number of components.');  end
  if size(A,dimT) ~= size(B,dimT), error('A and B have different number of test.'); end
  if size(A,dimT) ~= 1, error('not implemented yet... comming soon'); end

  if isvector( A ), A = A(:).'; end
  if isvector( B ), B = B(:).'; end

  if size(A,dimC) ~= 1, error('A and B must have 1 row.');  end
  
  TESTTYPE = 'unequal';
  [varargin,TESTTYPE] = parseargs(varargin,'UNEQual','$FORCE$',{'unequal',TESTTYPE});
  [varargin,TESTTYPE] = parseargs(varargin,  'EQual','$FORCE$',{  'equal',TESTTYPE});


  nA = size( A , dimS );
  nB = size( B , dimS );
  if nB == 0, TESTTYPE = [TESTTYPE,'1']; end
  
  sA2 = []; sA2_ = [];
  sB2 = []; sB2_ = [];
  switch lower(TESTTYPE)
    case   'equal1', t = T_equal1( A );
    case   'equal' , t = T_equal2( A , B );
    case 'unequal' , t = T_unequal2( A , B );
    otherwise, error('unknown or unable TEST TYPE.');
  end
  
  if nargout > 1
    PVALUETYPE = 'parametric';
    [varargin,PVALUETYPE] = parseargs(varargin,'PARAMetric','$FORCE$',{'parametric',PVALUETYPE});
    
    NP = 1e3;
    [varargin,i,NP] = parseargs(varargin,'PERMutation','$DEFS$',NP); if i, PVALUETYPE = 'permutation'; end
    
    switch [ lower(TESTTYPE) , '_' , lower(PVALUETYPE) ]
      case 'unequal_parametric',  [p,dfe] = welch_estimation( );
      case 'unequal_permutation'
        
        SPLITS = generaSPLITS( NP , nA , nB );

        p = NaN( size(t) );
        for tt = 1:numel(p)
          Pfcn = permutations( cat(dimS,A(:,:,tt),B(:,:,tt)) , SPLITS , @(varargin)abs(T_unequal2(varargin{:})) );
          p(1,1,tt) = Pfcn( abs( t(1,1,tt) ) );
        end
    end
    
  end
  
  function F = permutations( X , SPLITS , Tfcn )
    NP = size( SPLITS , 2 );
    ts = NaN( NP , 1 );
    for pp = 1:NP
      ts(pp) = feval( Tfcn , X(:,SPLITS(:,pp),1) , X(:,~SPLITS(:,pp),1) );
    end
    
    ts = sort( ts );
    ps = log10( linspace( 1 , (1/NP) , numel(ts) ).' );
    
    [kk,id1] = unique( ts , 'first' );
    [kk,id2] = unique( ts , 'last'  );
    e = ~( id1 == id2 ); id1 = id1(e); id2 = id2(e);
    for i = 1:numel(id1)
      ts( id2(i) ) = nextnum( ts(id2(i)) , 2 );
      ts( (id1(i)+1):(id2(i)-1) ) = NaN;
    end
    e = ~isnan( ts ); ts = ts(e); ps = ps(e);
    
    F = @(tvalue)pow10( Interp1D(ps,ts,tvalue,'linear','closest') );
  end
  
  function t = T_unequal2( A , B )
    sA2  = var(A,0,dimS);
    sA2_ = sA2 ./ nA;

    sB2  = var(B,0,dimS);
    sB2_ = sB2 ./ nB;

    t  = ( mean(A,dimS) - mean(B,dimS) ) ./ sqrt( sA2_ + sB2_ );
  end

  function [p,dfe] = welch_estimation()
    dfe = ( sA2_ + sB2_ ).^2 ./ ( sA2_.^2 ./ (nA-1) + sB2_.^2 ./ (nB-1) );
%     dfe = nA + nB - 2;
    p = 2 * tcdf( -abs(t) , dfe );
  end

end

