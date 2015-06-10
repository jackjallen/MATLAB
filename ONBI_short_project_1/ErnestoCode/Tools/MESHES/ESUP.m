function elements= ESUP( nodes , m )
% 
% elements = ESUP( nodes,mesh )
% 

if ~isfield( m,'ESUP' )
  m= CreateESUP( m );
end

elements = [];
for n=[nodes(:)]'
  elements= [elements  m.ESUP.el( m.ESUP.p(n)+1:m.ESUP.p(n+1)) ];
end
elements= unique( elements );
