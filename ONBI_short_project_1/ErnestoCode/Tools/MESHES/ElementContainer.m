function es= ElementCointainer( points , mesh )
%
% elements= ElementContainer( points , mesh )
%

if size(points,2) < 3
  points(:,3)=0;
end

es= zeros( size(points,1),1 );
ijk= SCAMCoord( points , mesh );

for p=1:size(points,1)
  elements= mesh.ELOC{ ijk(p,1),ijk(p,2),ijk(p,3) };
  for e=elements
    p1= mesh.xyz( mesh.tri(e,1),: );
    p2= mesh.xyz( mesh.tri(e,2),: );
    p3= mesh.xyz( mesh.tri(e,3),: );
    c1= sign( (p2(1)-p1(1))*(points(p,2)-p2(2)) - (p2(2)-p1(2))*(points(p,1)-p2(1)) );
    c2= sign( (p3(1)-p2(1))*(points(p,2)-p3(2)) - (p3(2)-p2(2))*(points(p,1)-p3(1)) );
    if c1*c2 < 0
      continue;
    end
    c3= sign( (p1(1)-p3(1))*(points(p,2)-p1(2)) - (p1(2)-p3(2))*(points(p,1)-p1(1)) );
    if ( c2 >= 0 & c3 >= 0 ) || ( c2 <= 0 & c3 <= 0 )
      es(p)= e;
      break;
    end
  end
end
