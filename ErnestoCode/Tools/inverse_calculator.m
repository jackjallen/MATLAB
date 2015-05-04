function [O,SSs] = inverse_calculator( x , level )

  if ~isa(x,'double')    , error( 'only a double real scalar number accepted.' ); end
  if ~isscalar( x )      , error( 'only a double real scalar number accepted.' ); end
  if imag( x ) ~= 0      , error( 'only a real number accepted.' ); end

  if nargin < 2, level = 100;      end
  
  
  FCN_families = {'basics','cos'};
  EPSs  = 2;
  
  invF = {}; strF = {}; dirF = {};

  if any( strcmp( 'basics' , FCN_families ) )

      invF{end+1} = @(x) x;           strF{end+1} = @(s) s;                               dirF{end+1} = @(x) x;
      invF{end+1} = @(x) 1./x;        strF{end+1} = @(s) ['1/',s];                        dirF{end+1} = @(x) 1./(x);
      invF{end+1} = @(x) -x;          strF{end+1} = @(s) ['-',s];                         dirF{end+1} = @(x) -(x);
      invF{end+1} = @(x) -x;          strF{end+1} = @(s) ['-(',s,')'];                    dirF{end+1} = @(x) -(x);

    for k = vect( [-1;1]* primes(5) )
      invF{end+1} = @(x)x.^(k);       strF{end+1} = @(s)[s,'^(1/',number2str(k,0),')'];   dirF{end+1} = @(x) x.^(1/k);
      invF{end+1} = @(x)x.^(k);       strF{end+1} = @(s)[s,'^1/',number2str(k,0)];        dirF{end+1} = @(x) x.^(1/k);
      invF{end+1} = @(x)x.^(k);       strF{end+1} = @(s)[s,'^(',number2str(1/k,0),')'];   dirF{end+1} = @(x) x.^(1/k);
      invF{end+1} = @(x)x.^(k);       strF{end+1} = @(s)[s,'^',number2str(1/k,0)];        dirF{end+1} = @(x) x.^(1/k);

      invF{end+1} = @(x)x.^(1/k);     strF{end+1} = @(s)[s,'^(',number2str(k,0),')'];     dirF{end+1} = @(x) x.^(k);
      invF{end+1} = @(x)x.^(1/k);     strF{end+1} = @(s)[s,'^',number2str(k,0)];          dirF{end+1} = @(x) x.^(k);
    end

  end


  if any( strcmp( 'basicsplus' , FCN_families ) )

    for k = vect( [-1;1]*setdiff( primes(20) , primes(5)  ) )
      invF{end+1} = @(x)x.^(k);       strF{end+1} = @(s)[s,'^(1/',number2str(k,0),')'];   dirF{end+1} = @(x) x.^(1/k);
      invF{end+1} = @(x)x.^(k);       strF{end+1} = @(s)[s,'^1/',number2str(k,0)];        dirF{end+1} = @(x) x.^(1/k);
      invF{end+1} = @(x)x.^(k);       strF{end+1} = @(s)[s,'^(',number2str(1/k,0),')'];   dirF{end+1} = @(x) x.^(1/k);
      invF{end+1} = @(x)x.^(k);       strF{end+1} = @(s)[s,'^',number2str(1/k,0)];        dirF{end+1} = @(x) x.^(1/k);

      invF{end+1} = @(x)x.^(1/k);     strF{end+1} = @(s)[s,'^(',number2str(k,0),')'];     dirF{end+1} = @(x) x.^(k);
      invF{end+1} = @(x)x.^(1/k);     strF{end+1} = @(s)[s,'^',number2str(k,0)];          dirF{end+1} = @(x) x.^(k);
    end

  end

  
  if any( strcmp( 'constants' , FCN_families ) )
    
    for k = unique( vect( [-1;1]*setdiff([ sqrt(2:5) , exp(-3:3) , 10.^(-5:5) ] , 1 ) ) )
      invF{end+1} = @(x)k*x;          strF{end+1} = @(s)[s,'/',number2str(k,0)];          dirF{end+1} = @(x) x/k;
      invF{end+1} = @(x)x/k;          strF{end+1} = @(s)[number2str(k,0),'*',s];          dirF{end+1} = @(x) k*x;
    end

  end
  
  if any( strcmp( 'sin' , FCN_families ) )
    
      invF{end+1} = @(x) sin(x);      strF{end+1} = @(s)['asin(' ,s,')'];                 dirF{end+1} = @(x) asin(x);
      invF{end+1} = @(x) asin(x);     strF{end+1} = @(s)['sin('  ,s,')'];                 dirF{end+1} = @(x) sin(x);
      invF{end+1} = @(x) sinh(x);     strF{end+1} = @(s)['asinh(',s,')'];                 dirF{end+1} = @(x) asinh(x);
      invF{end+1} = @(x) asinh(x);    strF{end+1} = @(s)['sinh(' ,s,')'];                 dirF{end+1} = @(x) sinh(x);
      invF{end+1} = @(x) sind(x);     strF{end+1} = @(s)['asind(',s,')'];                 dirF{end+1} = @(x) asind(x);
      invF{end+1} = @(x) asind(x);    strF{end+1} = @(s)['sind(',s,')'];                  dirF{end+1} = @(x) sind(x);
    for k = [-1 -1/2 1/2 1]
      invF{end+1} = @(x) asin(x)-k*pi;
      strF{end+1} = @(s)['sin(',s,'+',number2str(k*pi,25),')'];
      dirF{end+1} = @(x) sin(x+k*pi);

      invF{end+1} = @(x) asind(x)-k*180;
      strF{end+1} = @(s)['sind(',s,'+',number2str(k*180,25),')'];
      dirF{end+1} = @(x) sind(x+k*180);
    end

  end


  if any( strcmp( 'cos' , FCN_families ) )
    
      invF{end+1} = @(x) cos(x);      strF{end+1} = @(s)['acos(' ,s,')'];                 dirF{end+1} = @(x) acos(x);
      invF{end+1} = @(x) cos(x);      strF{end+1} = @(s)['-acos(',s,')'];                 dirF{end+1} = @(x) -acos(x);
      invF{end+1} = @(x) acos(x);     strF{end+1} = @(s)['cos('  ,s,')'];                 dirF{end+1} = @(x) cos(x);
      invF{end+1} = @(x) cosh(x);     strF{end+1} = @(s)['acosh(',s,')'];                 dirF{end+1} = @(x) acosh(x);
      invF{end+1} = @(x) cosh(-x);    strF{end+1} = @(s)['acosh(',s,')'];                 dirF{end+1} = @(x) acosh(x);
      invF{end+1} = @(x) acosh(x);    strF{end+1} = @(s)['cosh(' ,s,')'];                 dirF{end+1} = @(x) cosh(x);
      invF{end+1} = @(x) acosh(-x);   strF{end+1} = @(s)['-cosh(',s,')'];                 dirF{end+1} = @(x) -cosh(x);
      invF{end+1} = @(x) cosd(x);     strF{end+1} = @(s)['acosd(',s,')'];                 dirF{end+1} = @(x) acosd(x);
      invF{end+1} = @(x) cosd(x);     strF{end+1} = @(s)['-acosd(',s,')'];                dirF{end+1} = @(x) -acosd(x);
      invF{end+1} = @(x) acosd(x);    strF{end+1} = @(s)['cosd(',s,')'];                  dirF{end+1} = @(x) cosd(x);
    for k = [-2 -1.5 -1 -1/2 1/2 1 1.5 2]
      invF{end+1} = @(x) acos(x)-k*pi;
      strF{end+1} = @(s)['cos(',s,'+',number2str(k*pi,25),')'];
      dirF{end+1} = @(x) cos(x+k*pi);

      invF{end+1} = @(x) acosd(x)-k*180;
      strF{end+1} = @(s)['cosd(',s,'+',number2str(k*180,25),')'];
      dirF{end+1} = @(x) cosd(x+k*180);
    end

  end
  

  if any( strcmp( 'tan' , FCN_families ) )
    
      invF{end+1} = @(x) tan(x);      strF{end+1} = @(s)['atan(' ,s,')'];                 dirF{end+1} = @(x) atan(x);
      invF{end+1} = @(x) atan(x);     strF{end+1} = @(s)['tan('  ,s,')'];                 dirF{end+1} = @(x) tan(x);
      invF{end+1} = @(x) tanh(x);     strF{end+1} = @(s)['atanh(',s,')'];                 dirF{end+1} = @(x) atanh(x);
      invF{end+1} = @(x) atanh(x);    strF{end+1} = @(s)['tanh(' ,s,')'];                 dirF{end+1} = @(x) tanh(x);
      invF{end+1} = @(x) tand(x);     strF{end+1} = @(s)['atand(',s,')'];                 dirF{end+1} = @(x) atand(x);
      invF{end+1} = @(x) atand(x);    strF{end+1} = @(s)['tand(' ,s,')'];                 dirF{end+1} = @(x) tand(x);
    for k = [-2 -1.5 -1 -1/2 1/2 1 1.5 2]
      invF{end+1} = @(x) atan(x)-k*pi;
      strF{end+1} = @(s)['tan(',s,'+',number2str(k*pi,25),')'];
      dirF{end+1} = @(x) tan(x+k*pi);

      invF{end+1} = @(x) atand(x)-k*180;
      strF{end+1} = @(s)['tand(',s,'+',number2str(k*180,25),')'];
      dirF{end+1} = @(x) tand(x+k*180);
    end

  end


  if any( strcmp( 'cot' , FCN_families ) )
    
      invF{end+1} = @(x) cot(x);      strF{end+1} = @(s)['acot(' ,s,')'];                 dirF{end+1} = @(x) acot(x);
      invF{end+1} = @(x) acot(x);     strF{end+1} = @(s)['cot('  ,s,')'];                 dirF{end+1} = @(x) cot(x);
      invF{end+1} = @(x) coth(x);     strF{end+1} = @(s)['acoth(',s,')'];                 dirF{end+1} = @(x) acoth(x);
      invF{end+1} = @(x) acoth(x);    strF{end+1} = @(s)['coth(' ,s,')'];                 dirF{end+1} = @(x) coth(x);
      invF{end+1} = @(x) cotd(x);     strF{end+1} = @(s)['acotd(',s,')'];                 dirF{end+1} = @(x) acotd(x);
      invF{end+1} = @(x) acotd(x);    strF{end+1} = @(s)['cotd(' ,s,')'];                 dirF{end+1} = @(x) cotd(x);
    for k = [-2 -1.5 -1 -1/2 1/2 1 1.5 2]
      invF{end+1} = @(x) acot(x)-k*pi;
      strF{end+1} = @(s)['cot(',s,'+',number2str(k*pi,25),')'];
      dirF{end+1} = @(x) cot(x+k*pi);

      invF{end+1} = @(x) acotd(x)-k*180;
      strF{end+1} = @(s)['cotd(',s,'+',number2str(k*180,25),')'];
      dirF{end+1} = @(x) cotd(x+k*180);
    end

  end


  if any( strcmp( 'sec' , FCN_families ) )
    
      invF{end+1} = @(x) sec(x);      strF{end+1} = @(s)['asec(' ,s,')'];                 dirF{end+1} = @(x) asec(x);
      invF{end+1} = @(x) asec(x);     strF{end+1} = @(s)['sec('  ,s,')'];                 dirF{end+1} = @(x) sec(x);
      invF{end+1} = @(x) sech(x);     strF{end+1} = @(s)['asech(',s,')'];                 dirF{end+1} = @(x) asech(x);
      invF{end+1} = @(x) asech(x);    strF{end+1} = @(s)['sech(' ,s,')'];                 dirF{end+1} = @(x) sech(x);
      invF{end+1} = @(x) secd(x);     strF{end+1} = @(s)['asecd(',s,')'];                 dirF{end+1} = @(x) asecd(x);
      invF{end+1} = @(x) asecd(x);    strF{end+1} = @(s)['secd(' ,s,')'];                 dirF{end+1} = @(x) secd(x);
    for k = [-2 -1.5 -1 -1/2 1/2 1 1.5 2]
      invF{end+1} = @(x) asec(x)-k*pi;
      strF{end+1} = @(s)['sec(',s,'+',number2str(k*pi,25),')'];
      dirF{end+1} = @(x) sec(x+k*pi);

      invF{end+1} = @(x) asecd(x)-k*180;
      strF{end+1} = @(s)['secd(',s,'+',number2str(k*180,25),')'];
      dirF{end+1} = @(x) secd(x+k*180);
    end

  end


  if any( strcmp( 'csc' , FCN_families ) )
    
      invF{end+1} = @(x) csc(x);      strF{end+1} = @(s)['acsc(' ,s,')'];                 dirF{end+1} = @(x) acsc(x);
      invF{end+1} = @(x) acsc(x);     strF{end+1} = @(s)['csc('  ,s,')'];                 dirF{end+1} = @(x) csc(x);
      invF{end+1} = @(x) csch(x);     strF{end+1} = @(s)['acsch(',s,')'];                 dirF{end+1} = @(x) acsch(x);
      invF{end+1} = @(x) acsch(x);    strF{end+1} = @(s)['csch(' ,s,')'];                 dirF{end+1} = @(x) csch(x);
      invF{end+1} = @(x) cscd(x);     strF{end+1} = @(s)['acscd(',s,')'];                 dirF{end+1} = @(x) acscd(x);
      invF{end+1} = @(x) acscd(x);    strF{end+1} = @(s)['cscd(' ,s,')'];                 dirF{end+1} = @(x) cscd(x);
    for k = [-2 -1.5 -1 -1/2 1/2 1 1.5 2]
      invF{end+1} = @(x) acsc(x)-k*pi;
      strF{end+1} = @(s)['csc(',s,'+',number2str(k*pi,25),')'];
      dirF{end+1} = @(x) csc(x+k*pi);

      invF{end+1} = @(x) acscd(x)-k*180;
      strF{end+1} = @(s)['cscd(',s,'+',number2str(k*180,25),')'];
      dirF{end+1} = @(x) cscd(x+k*180);
    end

  end

  if any( strcmp( 'exp' , FCN_families ) )
    
      invF{end+1} = @(x) log(x);      strF{end+1} = @(s)['exp(',s,')'];                   dirF{end+1} = @(x) exp(x);
      invF{end+1} = @(x) exp(x);      strF{end+1} = @(s)['log(',s,')'];                   dirF{end+1} = @(x) log(x);

      invF{end+1} = @(x) log1p(x);    strF{end+1} = @(s)['expm1(',s,')'];                 dirF{end+1} = @(x) expm1(x);
      invF{end+1} = @(x) expm1(x);    strF{end+1} = @(s)['log1p(',s,')'];                 dirF{end+1} = @(x) log1p(x);

      invF{end+1} = @(x) log2(x);     strF{end+1} = @(s)['2^(',s,')'];                    dirF{end+1} = @(x) 2.^(x);
      invF{end+1} = @(x) log2(x);     strF{end+1} = @(s)['2^',s];                         dirF{end+1} = @(x) 2.^(x);
      invF{end+1} = @(x) 2.^(x);      strF{end+1} = @(s)['log2(',s,')'];                  dirF{end+1} = @(x) log2(x);

      invF{end+1} = @(x) log10(x);    strF{end+1} = @(s)['10^(',s,')'];                   dirF{end+1} = @(x) 10.^(x);
      invF{end+1} = @(x) log10(x);    strF{end+1} = @(s)['10^',s];                        dirF{end+1} = @(x) 10.^(x);
      invF{end+1} = @(x) 10.^(x);     strF{end+1} = @(s)['log10(',s,')'];                 dirF{end+1} = @(x) log10(x);

  end

  if any( strcmp( 'erf' , FCN_families ) )
    
      invF{end+1} = @(x) erf(x);      strF{end+1} = @(s)['erfinv(' ,s,')'];               dirF{end+1} = @(x) erfinv(x);
      invF{end+1} = @(x) erfinv(x);   strF{end+1} = @(s)['erf(' ,s,')'];                  dirF{end+1} = @(x) erf(x);

      invF{end+1} = @(x) erfc(x);     strF{end+1} = @(s)['erfcinv(' ,s,')'];              dirF{end+1} = @(x) erfcinv(x);
      invF{end+1} = @(x) erfcinv(x);  strF{end+1} = @(s)['erfc(' ,s,')'];                 dirF{end+1} = @(x) erfc(x);

  end
  
  
  
  S = double2str( x , 99 );
  SSs = { S }; SS = S;


  hWB = waitbar( 0 , S );
  for ff = 1:numel( invF )
    waitbar( ff/numel(invF) , hWB );
    try
      y = feval( invF{ff} , x );
      
      all_yy = unique( y + (-EPSs:EPSs)*eps(y) );
      all_yy = all_yy( feval( dirF{ff} , all_yy ) == x );
      for yy = all_yy
        SS = feval( strF{ff} , double2str(yy,99) );
        if ~any( strcmp(SSs,SS) )  &&  isequal( eval(SS) , x )
          if numel(SS) < numel(S), S = SS; end
          SSs = [ SSs ; SS ];
        end

        if level > 100
          SS = feval( strF{ff} , inverse_calculator( yy , level - 100 ) );
          if ~any( strcmp(SSs,SS) )  &&  isequal( eval(SS) , x )
            if numel(SS) < numel(S), S = SS; end
            SSs = [ SSs ; SS ];
          end
        end

      end
    end
  end
  close( hWB );

  
  
  SSs = uniquens( SSs );
  i = 0;
  while i < numel(SSs)
    i = i+1;
    SSo = SSs{i};

    SS = strrep(strrep(SSo,'(',''),')','');
    try, if ~any(strcmp(SS,SSs)) && isequal( eval( SS ) , x ), SSs = [ SSs ; SS ]; end; end

    SS = strrep(SSo,'--','');
    try, if ~any(strcmp(SS,SSs)) && isequal( eval( SS ) , x ), SSs = [ SSs ; SS ]; end; end

    SS = regexprep( SSo , '-(.*)/-(.*)', '$1/$2' );
    try, if ~any(strcmp(SS,SSs)) && isequal( eval( SS ) , x ), SSs = [ SSs ; SS ]; end; end

    SS = regexprep( SSo , '-(.*)\*-(.*)', '$1*$2' );
    try, if ~any(strcmp(SS,SSs)) && isequal( eval( SS ) , x ), SSs = [ SSs ; SS ]; end; end

    SS = regexprep( SSo , '-\(-(.*)\)' , '$1' );
    try, if ~any(strcmp(SS,SSs)) && isequal( eval( SS ) , x ), SSs = [ SSs ; SS ]; end; end

    SS = regexprep( SSo , '(.*)1/-(.*)', '$1-1/$2');
    try, if ~any(strcmp(SS,SSs)) && isequal( eval( SS ) , x ), SSs = [ SSs ; SS ]; end; end

    SS = regexprep( SSo , '(.*)\*1/(.*)', '$1/$2');
    try, if ~any(strcmp(SS,SSs)) && isequal( eval( SS ) , x ), SSs = [ SSs ; SS ]; end; end

    SS = regexprep( SSo , '(.*)/-(.*)', '-$1/$2');
    try, if ~any(strcmp(SS,SSs)) && isequal( eval( SS ) , x ), SSs = [ SSs ; SS ]; end; end

    SS = regexprep( SSo , '(.*)*-(.*)', '-$1*$2');
    try, if ~any(strcmp(SS,SSs)) && isequal( eval( SS ) , x ), SSs = [ SSs ; SS ]; end; end

    SS = regexprep( SSo , '(.*)\+-(.*)', '$1-$2');
    try, if ~any(strcmp(SS,SSs)) && isequal( eval( SS ) , x ), SSs = [ SSs ; SS ]; end; end

    SS = regexprep( SSo , '(.*)-\+(.*)', '$1-$2');
    try, if ~any(strcmp(SS,SSs)) && isequal( eval( SS ) , x ), SSs = [ SSs ; SS ]; end; end

    SS = regexprep( SSo , '(.*)\+0(.*)', '$1$2');
    try, if ~any(strcmp(SS,SSs)) && isequal( eval( SS ) , x ), SSs = [ SSs ; SS ]; end; end

    SS = regexprep( SSo , '(.*)-0(.*)', '$1$2');
    try, if ~any(strcmp(SS,SSs)) && isequal( eval( SS ) , x ), SSs = [ SSs ; SS ]; end; end

    SS = regexprep( SSo , '(.*)(\^1/)(.*)' , '$1/$3' );
    try, if ~any(strcmp(SS,SSs)) && isequal( eval( SS ) , x ), SSs = [ SSs ; SS ]; end; end

    try
    SS = regexp( SSo , '(.*\()(\d*\.?\d*)\+(\d*\.?\d*)(\).*)', 'tokens', 'once');
    SS = [ SS{1} number2str( eval(SS{2}) + eval(SS{3}) , 25 ) SS{4} ];
    try, if ~any(strcmp(SS,SSs)) && isequal( eval( SS ) , x ), SSs = [ SSs ; SS ]; end; end
    end
  
    try
    SS = regexp( SSo , '(.*\()(\d*\.?\d*)-(\d*\.?\d*)(\).*)', 'tokens', 'once');
    SS = [ SS{1} number2str( eval(SS{2}) - eval(SS{3}) , 25 ) SS{4} ];
    try, if ~any(strcmp(SS,SSs)) && isequal( eval( SS ) , x ), SSs = [ SSs ; SS ]; end; end
    end
  end
  SSs = uniquens( SSs );

  id = minInd( cellfun( 'prodofsize', SSs ) );
  if numel( SSs{id} ) < numel( S )
    S = SSs{id};
  end
  
  if nargout == 0
    disp(S)
  else
    O = S;
  end

end
