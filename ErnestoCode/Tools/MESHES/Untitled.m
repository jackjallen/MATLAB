
cd d:\work\MESHES\;
m= read_VTP('hipocampo.vtp');
m= CleanMesh( m );
m.triNORMALS= ComputeNormals( m );
bb=[-2.1116 -2.8479 -3.4339;1.2676 3.9105 3.3245];
[x,y,z]=ndgrid(linspace(bb(1,1),bb(2,1),256),linspace(bb(1,2),bb(2,2),128),linspace(bb(1,3),bb(2,3),256));

tic
  DT= GenerateImageDT( m , x , y , z );
toc

%%

% ts_vtk=[];
% for k=1:20
%   k
%   vtkClosestElement( m,points,k);
%   tic, vtkClosestElement( m,points,k); t= toc;
%   ts_vtk= [ts_vtk; k t];
% end

% ts=[];
for k=21:40
  k
  m= CreateLocators( m , k );
  ClosestElement( m,points);
  tic, ClosestElement( m,points); t= toc;
  ts= [ts; k t];
end


%%

ps= [];
for k=1:50
  m= CreateLocators( m , k );
  t= [];
  for j=1:100
    p= [rand(1,3)*5 ];
    tic; [ii,pp,dd]= ClosestPoint(m,p,k); t1= toc;
    tic; [d,i] = min( sqrt(sum( (m.xyz - repmat( p,size(m.xyz,1),1) ).^2 , 2))); t2= toc;

    t= [t t2/t1];
    if i ~= ii 
      ps= [ps;k p];
    end
  end
  disp(sprintf('%d --- %f', k,  exp( mean( log(t)))));
end

%%

ps = [];
t= [];
for j=1:1000
  p= randn(1,3)*20;
  tic; i= ClosestPoint(m,p); t1= toc;
  tic; [a,ii] = min( sum( (m.xyz - repmat( p,size(m.xyz,1),1) ).^2 , 2)); t2= toc;

  t= [t t2/t1];
  if i ~= ii
    ps= [ps;p];
  end
%   disp([j t2/t1]);
end
exp( mean( log(t)))

%%

p= randn(100,2); p(:,3)=0;

tic; i= ClosestPoint(m,p); toc

iii=[];
tic
for j=1:100
  [a,ii] = min( sum( (m.xyz - repmat( p(j,:),size(m.xyz,1),1) ).^2 , 2));
  iii(j,1)= ii;
end
toc



%%

for k= 1:20
  GenerateImageDT( m , x,y,z, .1 , k);
  k
  tic; DT= GenerateImageDT( m , x,y,z, .1 , k); toc
end

