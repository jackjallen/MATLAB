function mesh= MapMesh( mesh , type , varargin )
% 
% mesh_mapped= MapMeshSpherically( mesh , 'sph2cart'
%                                         'cart2sph'   
%                                         'pol2cart'   
%                                         'cart2pol'   
% 
% 'cart2sph'  [x y z] <-> [ theta phi r ]
% 'cart2pol'  [x y z] <-> [ theta z   r ]
%

if nargin > 2
  mesh= TransformMesh( mesh , varargin{:} );
end
  


if nargin > 1
  type= lower(type);
  fields= fieldnames( mesh );
  for f=1:size(fields,1)
    field= fields{f};
    if (strncmp( field, 'lmk',3) || strcmp( field, 'xyz') )
      if     ( strcmp( type,'cart2sph' ) )
        [t,p,r]= cart2sph( mesh.(field)(:,1) , mesh.(field)(:,2) , mesh.(field)(:,3) );
        mesh.(field)= [ t p r]; 
      elseif ( strcmp( type,'sph2cart' ) )
        [x,y,z]= sph2cart( mesh.(field)(:,1) , mesh.(field)(:,2) , mesh.(field)(:,3) );
        mesh.(field)= [ x y z];
      elseif ( strcmp( type,'cart2pol' ) )
        [t,r,z]= cart2pol( mesh.(field)(:,1) , mesh.(field)(:,3) , mesh.(field)(:,2) );
        mesh.(field)= [ t z r];
      elseif ( strcmp( type,'pol2cart' ) )
        [x,y,z]= cart2pol( mesh.(field)(:,1) , mesh.(field)(:,3) , mesh.(field)(:,2) );
        mesh.(field)= [ x y z];
      end
    end
  end
end
