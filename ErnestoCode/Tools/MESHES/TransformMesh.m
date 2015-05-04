function [M,H]= TransformMesh( M , varargin )
% [M,H]= TransformMesh( M , varargin )

  H = maketransform( varargin{:} );

  %now transform the points and landmarks fields
  if isfield( M , 'xyz' )
    M.xyz = transform( M.xyz , H , 'rows' );
  end
  if isfield( M , 'lmk' )
    M.lmk = transform( M.lmk , H , 'rows' );
  end
  
end





% % 
% % mesh= TransformMesh(  mesh   , 
% % 
% %                   , 'rotationUnits', [ default 'Degrees' , 'Radians' ]
% %                   , 'Center', [xc yc zc] or 'CenterMass' or 'BoundingBox'
% %                                          or 'Origin' (default)
% % 
% %                              , 'TranslateX' , tx,
% %                              , 'TranslateY' , ty,
% %                              , 'TranslateZ' , tz,
% %                              , 'RotateX' , rx
% %                              , 'RotateY' , ry
% %                              , 'RotateZ' , rz
% %                              , 'RelativeRotateZ' , rx
% %                       (rotate about Z axis about the center of the shape)
% %                              , 'RelativeRotateX' , ry
% %                              , 'RelativeRotateY' , rz
% %                              , 'Scale' , s
% %                              , 'RelativeScale' , s
% %                       (scale about the center of the shape)
% % 
% %                              , 'Pose3D', [ tx,ty,tz,rry,rrx,rrz,s ]
% %                              , 'TRANSFORM' ('H'), transform_matrix 
% % 
% %                              , 'CentertoCenterMass'
% %             translate the centermass of the mesh to 0,0,0
% %                              , 'CentertoBoundingBox'
% %             translate the center of the bounding box to 0,0,0
% %
% %                              , 'AxisAngle' , axis , angle
% %
% %                              , 'SimilarityAlign', from_points , to_points
% %
% %    to implement!
% %         align with respect of the moment axis
% %         align 3 points with others 3 points
% %         align 2 points with others 2 points
% %         rotation about axis angle representation
% %
% 
% 
% % xyz= mesh.xyz;
% 
% %some defaults
% units= 'degrees';
% center= 'origin';
% 
% i=1;
% while i <= length(varargin)
%   if     ( strcmp( lower(varargin{i}) , 'rotationunits' ) | strcmp( lower(varargin{i}),'u' ) )
%     switch lower(varargin{i+1})
%       case {'degrees','d'}
%         units= 'degrees';
%       case {'radians','r'}
%         units= 'radians';
%     end
%     i= i+2;
%   elseif ( strcmp( lower(varargin{i}) , 'center' ) | strcmp( lower(varargin{i}),'c' ) )
%     switch lower(varargin{i+1})
%       case {'centermass','cm'}
%         center= 'centermass';
%       case {'boundingbox','bb'}
%         center= 'boundingbox';
%       case {'origin','o'}
%         center= 'origin';
%       otherwise
%         center= varargin{i+1};
%     end
%     i= i+2;
%   else
%     i=i+1;
%   end
% end
% 
% if strcmp( 'radians' , units )
%   rotationscale= 1;
% elseif strcmp( 'degrees' , units )
%   rotationscale= pi/180;
% end
% 
% if strcmp( center , 'origin' )
%   center= [0 0 0];
% elseif strcmp( center , 'boundingbox' )
%   center= mean( [ min( mesh.xyz,[],1 ) ; max( mesh.xyz,[],1 ) ] );
% elseif strcmp(center,'centermass')
%   center= mean( mesh.xyz );
% end
% %center_transform
% c_tr= T( -center );
% %and its inverse, come back to the original position
% ic_tr= T( center );
% 
% 
% %initializing the transformation
% H= eye(4);
% 
% i=1;
% while i <= length(varargin)
%   switch lower( varargin{i} )
%     case {'translatex','tx'}
%       tx= varargin{i+1};
%       transform= T([tx,0,0]);
%       i=i+2;
%     case {'translatey','ty'}
%       ty= varargin{i+1};
%       transform= T([0,ty,0]);
%       i=i+2;
%     case {'translatez','tz'}
%       tz= varargin{i+1};
%       transform= T([0,0,tz]);
%       i=i+2;
%     case {'scale','s'}
%       s= varargin{i+1};
%       transform= diag([s s s 1]);
%       i=i+2;
%     case {'rotationx','rx'}
%       r= rotationscale*varargin{i+1};
%       transform= RX(r);
%       i=i+2;
%     case {'rotationy','ry'}
%       r= rotationscale*varargin{i+1};
%       transform= RY(r);
%       i=i+2;
%     case {'rotationz','rz'}
%       r= rotationscale*varargin{i+1};
%       transform= RZ(r);
%       i=i+2;
%     case {'relativescale','rs'}
%       s= varargin{i+1};
%       scale= diag([s s s 1]);
%       transform= ic_tr*scale*c_tr;
%       i=i+2;
%     case {'relativerotationx','rrx'}
%       r= rotationscale*varargin{i+1};
%       rotation= RX(r);
%       transform= ic_tr*rotation*c_tr;
%       i=i+2;
%     case {'relativerotationy','rry'}
%       r= rotationscale*varargin{i+1};
%       rotation= RY(r);
%       transform= ic_tr*rotation*c_tr;
%       i=i+2;
%     case {'relativerotationz','rrz'}
%       r= rotationscale*varargin{i+1};
%       rotation= RZ(r);
%       transform= ic_tr*rotation*c_tr;
%       i=i+2;
%     case {'pose3d','p3d'}
%       pose= varargin{i+1};
%       tx= pose(1);
%       ty= pose(2);
%       tz= pose(3);
%       ry= rotationscale*pose(4);
%       rx= rotationscale*pose(5);
%       rz= rotationscale*pose(6);
%       s = pose(7);
%       rotation= RZ(rz)*RX(rx)*RY(ry);
%       scale= diag([s s s 1]);
%       traslation= T([tx,ty,tz]);
%       transform= traslation*ic_tr*rotation*scale*c_tr;
%       i= i+2;
%     case {'transform','h'}
%       transform= varargin{i+1};
%       i=i+2;
%     case {'centertocentermass','ccm'}
%       centermass= mean( mesh.xyz );
%       transform= T( -centermass );
%       i=i+1;
%     case {'centertoboundingbox','cbb'}
%       centerbb=  mean( [ min( mesh.xyz,[],1 ) ; max( mesh.xyz,[],1 ) ] );
%       transform= T( -centerbb );
%       i=i+1;
%     case {'similarityalign','sa'}
%       from= varargin{i+1};
%       to  = varargin{i+2};
%       if     size( from , 1 ) == 1
%         transform= T(to)*T(-from);
%       elseif size( from , 1 ) == 2
%         Dto21   = to(2,:)-to(1,:);
%         Dfrom21 = from(2,:)-from(1,:);
%         %move from to origin
%         uno= T( -from(1,:) );
%         %scale
%         esc= norm( Dto21 )/ norm( Dfrom21 );
%         dos= diag([esc esc esc 1]);
%         %rotate
%         normal= cross( Dfrom21 , Dto21 );
%         normal= normal/norm(normal);
%         u= Dto21/norm(Dto21);                          %this part compute the angle
%         v= Dfrom21/norm(Dfrom21);                      %between two unitary vectors
%         angle= 2*atan2(norm(u-v), norm(u+v));      %----------------
%         
%         S=[  0          -normal(3)   normal(2) ;   %this part compute the rotation
%              normal(3)   0          -normal(1) ;   %of angle about the normal 
%             -normal(2)   normal(1)   0        ];   %----------------- 
%         rotation= eye(3) + sin(angle)*S + (1-cos(angle))*S*S;
%         tres= [ rotation zeros(3,1) ; 0 0 0 1];
%         %move to to
%         cuatro=T( to(1,:) );
% 
%         %concatenate the transformations
%         transform= cuatro*tres*dos*uno;
%       elseif size( from , 1 ) == 3
%         %the points 1 and 2 will be aligned and the third point
%         %will lie in the plane formed by the 3 points
%         Dto21   = to(2,:)-to(1,:);
%         Dfrom21 = from(2,:)-from(1,:);
%         %move from to origin
%         uno= T( -from(1,:) );
%         %scale
%         esc= norm( Dto21 )/ norm( Dfrom21 );
%         dos= diag([esc esc esc 1]);
%         %rotate1
%         normal= cross( Dfrom21 , Dto21 );
%         normal= normal/norm(normal);
%         u= Dto21/norm(Dto21);                          %this part compute the angle
%         v= Dfrom21/norm(Dfrom21);                      %between two unitary vectors
%         angle= 2*atan2(norm(u-v), norm(u+v));          %----------------
%         
%         S=[  0          -normal(3)   normal(2) ;       %this part compute the rotation
%              normal(3)   0          -normal(1) ;       %of angle about the normal 
%             -normal(2)   normal(1)   0        ];       %----------------- 
%         rotation= eye(3) + sin(angle)*S + (1-cos(angle))*S*S;
%         tres= [ rotation zeros(3,1) ; 0 0 0 1];
%         %now align the normal of both planes.
%         Dto31   = to(3,:)-to(1,:);
%         from= from';
%         from(4,:)=1;
%         from= tres*dos*uno*from;
%         from= (from(1:3,:)./repmat( from(4,:),3,1 ))';
%         Dfrom21 = from(2,:)-from(1,:);
%         Dfrom31 = from(3,:)-from(1,:);
%         Nto   = cross( Dto21, Dto31 );
%         Nto   = Nto/norm(Nto);
%         Nfrom = cross( Dfrom21, Dfrom31 );
%         Nfrom = Nfrom/norm(Nfrom);
%         
%         angle= 2*atan2(norm(Nto-Nfrom), norm(Nto+Nfrom));   %this part compute the angle
%                                                             %between two unitary vectors
%         Dto21= Dto21/norm(Dto21);
%         S=[  0          -Dto21(3)   Dto21(2) ;       %this part compute the rotation
%              Dto21(3)   0          -Dto21(1) ;       %of angle about the vector Dto21 
%             -Dto21(2)   Dto21(1)   0        ];       %----------------- 
%         rotation2= eye(3) + sin(angle)*S + (1-cos(angle))*S*S;
%         cuatro= [ rotation2 zeros(3,1) ; 0 0 0 1];
%         
%         %move to to
%         cinco=T( to(1,:) );
% 
%         %concatenate the transformations
%         transform= cinco*cuatro*tres*dos*uno;
%       elseif size( from , 1 ) > 3
%         %procrustes alignment of the points
%       end      
%       i=i+3;
%     case {'axisangle','aa'}
%       w= varargin{i+1};
%       w= w/norm(w);
%       r= rotationscale*varargin{i+2};
%       S=[  0    -w(3)  w(2) ;
%            w(3)  0    -w(1) ;
%           -w(2)  w(1)  0   ];
%       rotation= eye(3) + sin(r)*S + (1-cos(r))*S*S;
%       transform= [ rotation zeros(3,1) ; 0 0 0 1];
%       i=i+3;
%     otherwise
%       transform= eye(4);
%       i=i+1;
%   end
%   ic_tr(:,4) = transform*ic_tr(:,4);
%   c_tr = inv( ic_tr );
%   H= transform*H;
% end
% 
% %now transform the points and landmarks fields
% fields= fieldnames( mesh );
% for f=1:size(fields,1)
%   field= fields{f};
%   if (strncmp( field, 'lmk',3) || strcmp( field, 'xyz') )
%     %convert the points to 3d Homogeneous points
%     mesh.(field)= mesh.(field)';
%     mesh.(field)(4,:)=1;
%     %apply the transformation
%     mesh.(field)= H*mesh.(field);
%     %normalizing the homogeneous points and
%     %save the transformed points in the output as row-format
%     mesh.(field)= ( mesh.(field)(1:3,:)./repmat( mesh.(field)(4,:) , 3 ,1 ) )';
%   end
% end
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% %auxiliars functions
% function R= RX( r )
%   R=[ 1  0        0       0  ;
%       0  cos(r)  -sin(r)  0  ;
%       0  sin(r)   cos(r)  0  ;
%       0  0        0       1 ];
% function R= RY( r )
%   R=[  cos(r)  0  sin(r)  0  ;
%        0       1  0       0  ;
%       -sin(r)  0  cos(r)  0  ;
%        0       0  0       1 ];
% function R= RZ( r )
%   R=[ cos(r)  -sin(r)  0  0  ;
%       sin(r)   cos(r)  0  0  ;
%       0        0       1  0  ;
%       0        0       0  1 ];
% function translation= T( XYZ )
%   translation=[ 1  0  0  XYZ(1) ;
%                 0  1  0  XYZ(2) ;
%                 0  0  1  XYZ(3) ;
%                 0  0  0  1     ];
