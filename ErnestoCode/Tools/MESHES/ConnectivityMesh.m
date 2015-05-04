function c= ConnectivityMesh( m )
% 
% c= ConnectivityMesh( m )
% 

if ~isfield( m , 'PSUP'),  m= CreatePSUP(m); end
if ~isfield( m , 'ESUP'),  m= CreateESUP(m); end

region = 1;
nodes  = [1:size( m.xyz,1)]';

c = zeros(size(m.tri,1),1);
while numel(nodes)
  seed= nodes(1);
  npreg = 0;             %number of points in region
  newspreg  = [seed];        
  preg      = [seed];    %points in region
  while npreg ~= numel( preg )
    npreg = numel(preg);
    newspreg= setdiff( PSUP( newspreg,m ), preg);
    preg  = [preg newspreg];
  end
  nodes= setdiff( nodes,preg );
%   disp( numel( nodes ) );

%   tic
  preg = sort(preg);
  c( any(ismembc(m.tri,preg),2) ) = region;
%   toc
  
%   c( ismembc(preg,m.tri(:,1)) ) = region;
%   c( ismembc(preg,m.tri(:,2)) ) = region;
%   c( ismembc(preg,m.tri(:,3)) ) = region;
%   c( arrayfun( @(x) ismembc(x,preg) , m.tri(:,1) ) ) = region;
%   c( arrayfun( @(x) ismembc(x,preg) , m.tri(:,2) ) ) = region;
%   c( arrayfun( @(x) ismembc(x,preg) , m.tri(:,3) ) ) = region;
%   c( ESUP(preg,m),1 )= region;

  region= region+1;
end
