function tmp_name = tmpname( pre , n )
%TMPNAME Get temporary file.
%   TEMPNAME returns a unique name, starting with the directory returned by 
%   TEMPDIR, suitable for use as a temporary file.
%
%   See also TEMPDIR.

  if nargin < 1
    pre = 8;
  end
  if nargin < 2
    if isnumeric( pre )
      n = pre;
      pre = 'temp_';
    else
      n = 8;
    end
  end
  
  
  pre = strrep( strrep( pre , '\' , filesep ) , '/' , filesep );
  
  
  lastdot = find( pre == '.' , 1 , 'last');
  if isempty( lastdot )
    ext = '';
  else
    ext = pre(lastdot:end);
    pre = pre(1:lastdot-1);
  end
  if strcmp(ext,'.')
    ext = '';
  end
  
  asdir = false;
  if pre(end) == filesep
    asdir = true;
    pre = pre(1:end-1);
  end
  

  temporal_dir = tempdir;
  if isunix
    [status,result] = system( ['ls -a /scratch'] );
    if ~status
      temporal_dir = fullfile( '/scratch/' , getenv('USER') );
    end
  end
  
  try,
  if isequal( pre(1:numel(temporal_dir)) , temporal_dir )
    pre = pre(numel(temporal_dir)+1:end);
  end
  end
  
  it = 0;
  while 1
    if n > 0
      uu = [];
      if isempty(uu), uu = strrep( uuid , '-' , '_' );                      end
      if isempty(uu), uu = strrep(char(java.util.UUID.randomUUID),'-','_'); end
      if isempty(uu), uu = num2str(feature('timing','cpucount'));           end

      if isnumeric( n )  &&  isscalar( n )
        uu = uu(1:n);
      elseif ischar( n )
        uu = uu(1:numel(n));
        m = isstrprop( uu , 'digit' );
        uu( m ) = uu( m ) - '1' + 'b';
      end
    elseif n < 0
      uu = sprintf('%0*d',-n,it);
    elseif n==0
      uu = sprintf('%d',it);
    end
      
    tmp_name = fullfile( temporal_dir , [ pre , uu , ext ] );

    [status,result] = system( ['ls -a ' tmp_name ] );
    if status
      break;
    end
    it = it+1;
    if ( n>0 && it > 100 ) || ( n==0 && it > 1e6 ) || ( n<0 && it > 10^(-n)-1 )
      error('esta lleno tempdir de files??');
    end
  end
  
  if asdir
    tmp_name = [tmp_name filesep];
  end
  
  tmp_name = strrep( tmp_name , filesep , '\' );
  tmp_name = regexprep( tmp_name , '\\+' , '\' );
  tmp_name = strrep( tmp_name , '\' , filesep );
  
  if asdir
    temporal_dir = tmp_name;
  else
    temporal_dir = fileparts( tmp_name );
  end
  if ~isdir( temporal_dir ), mkdir( temporal_dir ); end
  
end
