function x = trand( varargin )
% 
% trand - true rand function
% 

  %just to know the size!
  x = rand( varargin{:} );
  
  N = 10000;
  for i = 0:floor(numel(x)/N)-1
    x( i*N + (1:N) ) = sscanf( ...
      urlread(sprintf('https://www.random.org/integers/?num=%d&min=-1000000000&max=1000000000&col=1&base=10&format=plain&rnd=new',N)) ,...
      '%g' );
  end
  
  if numel(x) < N, i = -1; end
  x( ( (i+1)*N+1 ):numel(x) ) = sscanf( ...
      urlread(sprintf('https://www.random.org/integers/?num=%d&min=-1000000000&max=1000000000&col=1&base=10&format=plain&rnd=new',rem(numel(x),N))) ,...
                          '%g' );
  x = ( x + 1000000000 )/2000000000;
  
% urlread('http://www.random.org/sequences/?min=1&max=10&col=1&format=plain
% &rnd=new')
% 
% urlread('http://www.random.org/keno/?tickets=1&numbers=20')
% 
% urlread('http://www.random.org/strings/?num=1&len=20&digits=off&upperalpha=on&loweralpha=off&unique=off&format=plain&rnd=new')
% 
% http://www.random.org/cgi-bin/randbitmap?format=png&width=200&height=200&zoom=1

end
