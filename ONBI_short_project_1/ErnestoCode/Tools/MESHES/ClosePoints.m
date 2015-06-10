function points= ClosePoints( pto , mesh )
% 
% points= ClosePoints( pto , mesh )
% 

  D= min([ mesh.SCAM.Dx mesh.SCAM.Dy mesh.SCAM.Dz]);
  N= 0;
  points= [];
  while isempty( points )
    ijk1= SCAMCoord( pto - N*D , mesh );
    ijk2= SCAMCoord( pto + N*D , mesh );
    for i=ijk1(1):ijk2(1)
      for j=ijk1(2):ijk2(2)
        for k=ijk1(3):ijk2(3)
          points= [ points mesh.PLOC{i,j,k} ];
        end
      end
    end
    N=N+1;
  end
  
  points= unique(points);
  points= points(:)';
