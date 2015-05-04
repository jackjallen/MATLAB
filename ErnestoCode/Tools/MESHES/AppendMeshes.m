function m1= AppendMeshes( m1 , m2 , varargin )
% 
% m= AppendMeshes( m1 , m2 )
% 

if ~isfield(m1,'tri'), m1=struct('xyz',[],'tri',[]);end;
if ~isfield(m2,'tri'), m2=struct('xyz',[],'tri',[]);end;

if isempty(m1) 
  m1 = m2;
  return;
end
if isempty(m2)
  return;
end


np1= size( m1.xyz,1);
nt1= size( m1.tri,1);
np2= size( m2.xyz,1);
nt2= size( m2.tri,1);

fields= unique([fieldnames( m1 ) ; fieldnames( m2 )]);
for f=1:size(fields,1)
  field= fields{f};

  if strcmp( field , 'xyz')
    m1.xyz= [ m1.xyz ; m2.xyz ];
    continue;
  end

  if strcmp( field , 'tri')
    m1.tri= [ m1.tri ; m2.tri + np1 ];
    continue;
  end

  if ( strncmp( field , 'xyz' , 3 ) )
    if isfield(m2,field)
      m1.(field)(np1+1:np1+np2,:)= m2.(field);
    else
      m1.(field)(np1+np2,1)= 0;
    end
    continue;
  end

  if ( strncmp( field , 'tri' , 3 ) )
    if isfield(m2,field)
      m1.(field)(nt1+1:nt1+nt2,:)= m2.(field);
    else
      m1.(field)(nt1+nt2,1)= 0;
    end
    continue;
  end

 
  if      ( isfield(m1,field) & isfield(m2,field) )
    m1.(field)= [ m1.(field) ; m2.(field) ];
  elseif  ( ~isfield(m1,field) & isfield(m2,field) )
    m1.(field)= m2.(field);
  end

end


