function [Mout,Min,Nout,Nin]=Intersection3D(M,N);
% [Mout,Min,Nout,Nin]=Intersecton3D(M,N);

M.trilabel=zeros(size(M.tri,1),1);
N.trilabel=ones(size(M.tri,1),1);

difMN=vtkBooleanOperationPolyDataFilter(M,N,'SetOperation',2,'ReorientDifferenceCellsOff',[]);
difNM=vtkBooleanOperationPolyDataFilter(N,M,'SetOperation',2,'ReorientDifferenceCellsOff',[]);

Mout=difMN; Nin=difMN;
Min=difNM; Nout=difNM;

fields= fieldnames( M ) ; 
for f=1:size(fields,1)
  field= fields{f};
  if strncmp(field,'tri',3) && ~strcmp(field,'trilabel')
      Mout.(field)(Mout.trilabel==0,:)=[];
      Min.(field)(Min.trilabel==0,:)=[];
  end;
end;

fields= fieldnames( N ) ;
for f=1:size(fields,1)
    field= fields{f}; 
    if strncmp(field,'tri',3) && ~strcmp(field,'trilabel')
      Nout.(field)(Nout.trilabel==1,:)=[];
      Nin.(field)(Nin.trilabel==1,:)=[];
  end;
end;

Mout=rmfield(Mout,'trilabel');
Min=rmfield(Min,'trilabel');
Nout=rmfield(Nout,'trilabel');
Nin=rmfield(Nin,'trilabel');