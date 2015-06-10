function elements= CloseElements( pto , mesh )
% 
% elements= CloseElements( pto , mesh )
% 


  D= min([ mesh.SCAM.Dx mesh.SCAM.Dy mesh.SCAM.Dz]);
  N= 0;

  elements= [];
  while isempty( elements )
    ijk1= SCAMCoord( pto - N*D , mesh );
    ijk2= SCAMCoord( pto + N*D , mesh );
    for i=ijk1(1):ijk2(1)
      for j=ijk1(2):ijk2(2)
        for k=ijk1(3):ijk2(3)
          elements= [ elements mesh.ELOC{i,j,k} ];
        end
      end
    end
    N=N+1;
  end

  for i=ijk1(1):ijk2(1)
    for j=ijk1(2):ijk2(2)
      for k=ijk1(3):ijk2(3)
        elements= [ elements mesh.ELOC{i,j,k} ];
      end
    end
  end
  
  elements= unique(elements);
  elements= elements(:)';
