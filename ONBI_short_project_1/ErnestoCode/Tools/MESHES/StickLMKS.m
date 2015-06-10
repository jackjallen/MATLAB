function m= StickLMKS( m )
% 
% m= StickLMKS( m )
% 

m= CreateELOC( CreateSCAM( m,1) );

fields= fieldnames( m );
for f=1:size(fields,1)
  field= fields{f};
  if ( strncmp( field, 'lmk',3) )
    [e,lmks]= ClosestElement( m.(field) , m );
    m.(field)= lmks;
  end
end
