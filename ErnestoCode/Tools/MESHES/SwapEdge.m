function mesh= SwapRandomEdge( mesh , node1 , node2 )
% 
% mesh= SwapRandomEdge( mesh , node1 , node2 )
% 

  if (~isfield( mesh,'ESUP') )
    mesh= CreateESUP( mesh );
  end
  if (~isfield( mesh,'PSUP') )
    mesh= CreatePSUP( mesh );
  end

  elements= intersect( ESUP( node1,mesh) , ESUP( node2,mesh ) );
  if( numel(elements)== 2 )
    e1= elements(1);
    e2= elements(2);
    %node3 belong to e1
    node3= setdiff( mesh.tri(e1,:),[node1 node2] );
    %node4 belong to e2
    node4= setdiff( mesh.tri(e2,:),[node1 node2] );

    mesh.tri( e1,: )= [ node3 node4 node1 ];
    mesh.tri( e2,: )= [ node4 node3 node2 ];

    %update the PSUP structure
    mesh= remove_point( mesh , node1 , node2);
    mesh= remove_point( mesh , node2 , node1);
    mesh= add_point( mesh , node3 , node4);
    mesh= add_point( mesh , node4 , node3);

    %update the ESUP structure
    mesh= remove_element( mesh , node1 , e1);
    mesh= remove_element( mesh , node1 , e2);
    mesh= remove_element( mesh , node2 , e1);
    mesh= remove_element( mesh , node2 , e2);
    mesh= add_element( mesh , node1 , e1);
    mesh= add_element( mesh , node2 , e2);
    mesh= add_element( mesh , node3 , e2);
    mesh= add_element( mesh , node4 , e1);

  else
    fprintf(' -- imposible to swap\n');
  end
end


function mesh= remove_point ( mesh , p1 , p2 )
  mesh.PSUP.points= [ mesh.PSUP.points( 1:mesh.PSUP.p(p1)) ...
    setdiff( mesh.PSUP.points( mesh.PSUP.p(p1)+1:mesh.PSUP.p(p1+1) ), p2 ) ...
    mesh.PSUP.points( mesh.PSUP.p(p1+1)+1:end ) ];
  mesh.PSUP.p(p1+1:end) = mesh.PSUP.p(p1+1:end)-1;
end
function mesh= add_point ( mesh , p1 , p2 )
  mesh.PSUP.points= [ mesh.PSUP.points( 1:mesh.PSUP.p(p1)) ...
    sort([ mesh.PSUP.points( mesh.PSUP.p(p1)+1:mesh.PSUP.p(p1+1) ) p2] ) ...
    mesh.PSUP.points( mesh.PSUP.p(p1+1)+1:end ) ];
  mesh.PSUP.p(p1+1:end) = mesh.PSUP.p(p1+1:end)+1;
end

function mesh= remove_element( mesh , p , e )
  mesh.ESUP.el= [ mesh.ESUP.el( 1:mesh.ESUP.p(p) ) ...
    setdiff( mesh.ESUP.el( mesh.ESUP.p(p)+1:mesh.ESUP.p(p+1) ), e ) ...
    mesh.ESUP.el( mesh.ESUP.p(p+1)+1:end ) ];
  mesh.ESUP.p(p+1:end) = mesh.ESUP.p(p+1:end)-1;
end
function mesh= add_element( mesh , p , e )
  mesh.ESUP.el= [ mesh.ESUP.el( 1:mesh.ESUP.p(p)) ...
    sort([ mesh.ESUP.el( mesh.ESUP.p(p)+1:mesh.ESUP.p(p+1) ) e] ) ...
    mesh.ESUP.el( mesh.ESUP.p(p+1)+1:end ) ];
  mesh.ESUP.p(p+1:end) = mesh.ESUP.p(p+1:end)+1;
end


