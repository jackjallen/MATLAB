function mesh= SwapRandomEdge( mesh , node1 , node2 )
% 
% mesh= SwapRandomEdge( mesh )
% 

mesh= CreateESUP( mesh );
mesh= CreatePSUP( mesh );

if nargin < 2
  node1= ceil( rand(1)*size(mesh.xyz,1) );
end
fprintf('node1: %6d ' , node1 );

if nargin < 3
  psunode1= PSUP( node1 , mesh );
  node2= psunode1( ceil( rand(1)*numel( psunode1 ) ) );
end
fprintf(' --  node2: %6d ' , node2 );

elements= intersect( ESUP( node1,mesh) , ESUP( node2,mesh ) );
if( numel(elements)== 2 )
  fprintf(' -- elements: %6d %6d ', elements );
  nodes34= mesh.tri( elements,: );
  nodes34= setdiff( nodes34(:) , [node1 node2] );
  p1= mesh.xyz( node1,: );
  p2= mesh.xyz( node2,: );
  p3= mesh.xyz( nodes34(1),: );
  p4= mesh.xyz( nodes34(2),: );
  
  
  if intersectlines(p1,p2,p3,p4)
    mesh.tri( elements(1),: )= [ nodes34(1) nodes34(2) node1 ];
    mesh.tri( elements(2),: )= [ nodes34(1) nodes34(2) node2 ];
    mesh= CreateESUP( mesh );
    mesh= CreatePSUP( mesh );
    fprintf(' -- OK\n');
  else
    fprintf(' -- imposible to swap\n');
  end
else
  fprintf(' -- imposible to swap\n');
end


function i= intersectlines( p1 , p2,p3,p4)
  A=(p3(1)-p1(1))/(p2(1)-p1(1));
  B=(p4(1)-p3(1))/(p2(1)-p1(1));
  C=(p1(2)-p3(2))/(p4(2)-p3(2));
  D=(p2(2)-p1(2))/(p4(2)-p3(2));

  ta= (A+B*C)/(1-B*D);
  tb= C+D*ta;
  if( ta < 1 && ta > 0 && tb < 1 && tb > 0)
    i=1;
  else
    i=0;
  end
  
