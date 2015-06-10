function points= PSUP( nodes , m )
% 
% points = ESUP( nodes,mesh )
% 

if ~isfield( m,'PSUP' )
  m= CreatePSUP( m );
end

points = [];
for n=[nodes(:)]'
  points= [points  m.PSUP.points( m.PSUP.p(n)+1:m.PSUP.p(n+1)) ];
end
% points= unique( setdiff(points,nodes) );
points= unique( points );
