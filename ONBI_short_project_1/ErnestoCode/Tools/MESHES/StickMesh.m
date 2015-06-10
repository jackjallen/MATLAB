function original= StickMesh( original , final )
% 
% mesh= StickMesh( original , final )
% 

final= CreateELOC( CreateSCAM( final, 1) );
[ e, new_points ]= ClosestElement( original.xyz , final );

original.xyz= new_points;

