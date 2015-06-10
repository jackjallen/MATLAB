function varargout = parseargs( args , varargin )
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
%
%{

args = {'pl','v',4,-5,'method',pi+4};
[args,i,value1,value2]  = parseargs( args , 'ValueUnits' ,'$CHECK$',@(x) x>0 ,'$DEFS$',0,1)   -> no encuentra 'vu' entonces devuelve los defaults

[args,i,value1,value2]  = parseargs( args , 'ValueUnits',{'valueunits','methods','plot'},'$CHECK$',@(x) x>0 ,'$DEFS$',0,1)  -> encuentra 'v', que como es nonambiguous para las palabras que se difinieron, lo utiliza como valueunits

[args,i,value1,value2]  = parseargs( args , 'ValueUnits' ,'$CHECK$',@(x) x>0 ,'$DEFS$',0,1)
[args,hayMeth,meth]     = parseargs( args , 'Method' , '$FORCE$')
[args,pl,plidx]               = parseargs( args , 'plot',{'valueunits','methods','plot'},'$FORCE$',{1,'noPlot'})
a                       = parseargs( args , 'plot',{'valueunits','methods','plot'},'$FORCE$',{1,'noPlot'})
a                       = parseargs( args , 'plot',{'valueunits','methods','plot'},'$FORCE$',1)
%}

  try,

  if nargout > 1
    nvalues = nargout - 2;
  else
    nvalues = 1;
  end

  defs_id  = find( cellfun( @(s) strcmp(s,'$DEFS$' ) , varargin ) );
  if ~isempty( defs_id )
    if ~isscalar( defs_id ), error('muchos $DEFS$ .'); end
    if defs_id + nvalues ~= numel( varargin )
      error('after ''$DEFS$'' a default value/values are expected and it has to be at the end!!!.');
    end
    values   = varargin( defs_id + (1:(nargout-2)) );
    varargin = varargin( 1:(defs_id-1) );
  else
    values = cell(1,nvalues);
  end

  
  force_id = find( cellfun( @(s) strcmp(s,'$FORCE$') , varargin ) );
  if ~isempty( force_id )
    if ~isscalar( force_id ), error('muchos $FORCE$s .'); end
    if force_id + 1 > numel( varargin )
      %error('after ''$FORCE$'' a return or a cell value is expected.');
      force  = { true , false };
    else
      force  = varargin{ force_id + 1};
    end
    varargin = varargin( [ 1:(force_id-1)   (force_id+2):end ] );
  else
    force = [];
  end
  if iscell( force )
    if numel( force ) > 2
      error('si force is cell has to be have 2 elements!!');
    elseif numel( force ) == 1
      force{2} = false;
    elseif numel( force ) == 0
      force = {true,false};
    end
  elseif ~isempty( force )
    force = {force,false};
  end


  check_id = find( cellfun( @(s) strcmp(s,'$CHECK$') , varargin ) );
  if ~isempty( check_id )
    if ~isscalar( check_id ), error('muchos $CHECK$ .'); end
    if check_id + 1 > numel( varargin )
      error('after ''$CHECK$'' a checking is expected.');
    end
    check  = varargin{ check_id + 1};
    varargin = varargin( [ 1:(check_id-1)   (check_id+2):end ] );
  else
    check = [];
  end
  
  
  words = {};
  words_id = find( cellfun( @(s) iscell(s) , varargin ) );
  if ~isempty( words_id )
    if ~isscalar( words_id ), error('muchos  VALIDS_WORDS .'); end
    words = varargin{ words_id };
    varargin = varargin( [ 1:(words_id-1)   (words_id+1):end ] );
  end
  
  
  if any( cellfun( @(s) ~isempty( regexp( s , '^\$.*' , 'ONCE' ) ) , varargin( cellfun( @(s) ischar(s) , varargin ) ) ) )
    error( 'invalid $option' );
  end
  
  
  keys = varargin;
  found = 0;
  a_id = 1;
  while a_id <= numel(args)
    a = args{a_id};
    if ~ischar(a), a_id = a_id + 1; continue; end
    try, a = validatestring( a , words ); end
    
    for k_id = 1:numel(keys)
      k = keys{k_id};

      if ~strcmpi( a , k )  &&  ~isequal( upper(a)  ,  k(isstrprop(k,'upper')) )
        continue; 
      end
      
      found  = a_id;
      values = args( a_id + (1:nvalues) );
      
      if ~isempty( check ),
        try
          if any( ~cellfun( @(v) check(v) , values ) )
            fprintf('invalid arguments , values do not pass checkFcn\n');
            error('invalid arguments , values do not pass checkFcn');
          end
        catch
          error('error checking arguments');
        end
      end
      
      args   = args( [ 1:(a_id-1)  ,  (a_id+nvalues+1):end ] );
      a_id = a_id - 1;
      break;
    end
    a_id = a_id + 1;
%     if found, break; end
  end
  
  
  if ~isempty( force )
    if found
      if iscell( force )
        found = force{1};
      else
        found = force;
      end
    else
      if iscell( force )
        found = force{2};
      else
        found = false;
      end
    end    
  end


  if nargout > 1
    varargout{ 1 } = args;
    varargout{ 2 } = found;
    varargout( 2 + (1:nvalues) ) = values;
  else
    varargout{1} = found;
  end
  
  
  catch LE
    
    disperror(LE)
    
    
    error('parseargs ha cambiado un poco.... hay una version parseargs_old');
    
  end
  
end
