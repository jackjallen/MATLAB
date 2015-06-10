function bb = BBMesh( m )
% 
% bb = BBMesh( m )
% 

bb= [ min(m.xyz,[],1) ; max(m.xyz,[],1) ];
