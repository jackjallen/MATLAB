f= read_OBJ( 'd:\pru\AV0SP_001_3D1_Frontal.obj');
s= 1:10:100;
tiempos= [];
for ss= s
  tic
  CleanMesh( f , 'scam', ss )
  t= toc
  tiempos= [ tiempos t ];
end

