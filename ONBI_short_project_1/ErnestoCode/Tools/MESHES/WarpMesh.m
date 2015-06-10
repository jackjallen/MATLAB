function mesh= WarpMesh( mesh , o , t , varargin)
% 
% mesh= WarpMesh( mesh , o , t , 
%                            'NumberSpatialDimensions', [2 3]
%                            'Plot' 
%                            'BASIS', [ 'R','RLOGR' ]
% 

zz= parseargs( varargin,'plot','p' );

fields= fieldnames( mesh );
for f=1:size(fields,1)
  field= fields{f};
  if (strncmp( field, 'lmk',3) || strcmp( field, 'xyz') )
    mesh.(field)= WarpPoints( mesh.(field) , o , t , varargin{[1:zz-1 zz+1:end]} );
  end
end

if zz
  WarpPoints( o , o , t , varargin{1:end} );
end
