function x = trandn( varargin )
% 
% trandn - true randn function
% 

  x = rand( varargin{:} )*0;
  
  N = 10000;
  for i = 0:floor(numel(x)/N)-1
    x( i*N + (1:N) ) = sscanf( ...
      urlread(sprintf('https://www.random.org/gaussian-distributions/?num=%d&mean=0.0&stdev=1.0&dec=20&col=1&notation=scientific&format=plain&rnd=new',N)) ,...
      '%g' );
  end
  
  if numel(x) < N, i = -1; end
  x( ( (i+1)*N+1 ):numel(x) ) = sscanf( ...
      urlread(sprintf('https://www.random.org/gaussian-distributions/?num=%d&mean=0.0&stdev=1.0&dec=20&col=1&notation=scientific&format=plain&rnd=new',rem(numel(x),N))) ,...
                          '%g' );
  

end

