function plotMesh( mesh , varargin )
%
% plotmesh( mesh , 'Text'
%                  'NumberPoints'
%                  'NumberElements'
%                  'MarkerPoint'   , Symbol to mark the points
%                  'TriangleColor' , [r g b] or 'none'   ('FaceColor')              
%                  'EdgeColor'     , [r g b] or 'none'                 
%                  'PointData'     , data_values_on_points
%                  'TriangleData'  , data_values_on_triangles  ('FaceData')
%
%



if isfield( mesh , 'vertices' ) && isfield( mesh , 'faces' )
  mesh = struct( 'xyz',mesh.vertices , 'tri', mesh.faces );
end



lmkid=1;
lmkls= [ '.b' ; '.g' ; '.r' ; '.c' ; '.m' ; '.y' ; '.k' ];

markersymbol= 'none';

elements= 1:size(mesh.tri,1);
texto= 'no';
numberpoints= 'no';
numberelements= 'no';
facecolor= [.6 .75 .75];
edgecolor= [ 0 0 0 ];
pointdata= [];
triangledata= [];
                 
[varargin,zz, elements      ] = parseargs( varargin,'elements'                   ,'$DEFS$', elements     );          
[varargin,zz, facecolor     ] = parseargs( varargin,'TriangleColor' ,'FaceColor' ,'$DEFS$', facecolor    );  
[varargin,zz, edgecolor     ] = parseargs( varargin,'EdgeColor'                  ,'$DEFS$', edgecolor    );      
[varargin,zz, triangledata  ] = parseargs( varargin,'TriangleData','FaceData'    ,'$DEFS$', triangledata );
[varargin,zz, pointdata     ] = parseargs( varargin,'PointData'                  ,'$DEFS$', pointdata    );
[varargin,zz, markersymbol  ] = parseargs( varargin,'MarkerPoint'                ,'$DEFS$', markersymbol );

[varargin, texto            ] = parseargs( varargin,'text'                       ,'$FORCE$',{'yes',texto} );
[varargin, numberpoints     ] = parseargs( varargin,'numberpoints'               ,'$FORCE$',{'yes',numberpoints} );
[varargin, numberelements   ] = parseargs( varargin,'numberelements'             ,'$FORCE$',{'yes',numberelements} );

    
if strcmp( texto,'yes' )
  numberpoints= 'yes';
  numberelements= 'yes';
end


% fix the data if any triangle is degenerate...
if ~isempty(mesh.tri)
    mesh= FixMesh( mesh );
end;

rotate3d off;
newplot;
bloquea= ishold;
if ~bloquea, cla; end

if      size(pointdata,2)== 1
  patch( 'vertices', mesh.xyz , 'faces', mesh.tri(elements,:), ...
        'facecolor', 'interp' , ...
        'edgecolor', edgecolor , ...
        'facevertexcdata', pointdata ,varargin{:});
  patch( 'vertices', mesh.xyz , 'faces', vec(1:size(mesh.xyz,1)), ...
        'facecolor', 'interp' , ...
        'Marker' , markersymbol , ...
        'MarkerFaceColor', 'flat', ...
        'Cdata', pointdata, ...
        varargin{:});
    
elseif  size(triangledata,2)== 1
  if ( islogical(elements) && size( triangledata , 1 ) ~= nnz( elements ) ) || ( size( triangledata , 1 ) ~= numel( elements ) )
    triangledata = triangledata( elements );
  end
  
  patch( 'vertices', mesh.xyz , 'faces', mesh.tri(elements,:), ...
        'facecolor', 'flat' , ...
        'Marker' , markersymbol , ...
        'edgecolor', edgecolor , ...
        'facevertexcdata', triangledata ,varargin{:});
    
else
  patch( 'vertices', mesh.xyz , 'faces', mesh.tri(elements,:), ...
         'facecolor', facecolor , ...
         'edgecolor', edgecolor ,varargin{:});
    hold on;
  if ~strcmp(markersymbol,'none')  
      patch( 'vertices', mesh.xyz , 'faces', vec(1:size(mesh.xyz,1)), ...
        'Marker' , markersymbol , varargin{:});
  end;
end

hold on;
axis image;
xlabel('X');
ylabel('Y');
zlabel('Z');

if      size(pointdata,2)== 3
  quiver3(mesh.xyz(:,1),mesh.xyz(:,2),mesh.xyz(:,3), ...
        pointdata(:,1),pointdata(:,2),pointdata(:,3),'b');
end
if      size(triangledata,2)== 3
  centers=[ mean( reshape(mesh.xyz(mesh.tri(:),1),size(mesh.tri)),2 ) ...
            mean( reshape(mesh.xyz(mesh.tri(:),2),size(mesh.tri)),2 ) ...
            mean( reshape(mesh.xyz(mesh.tri(:),3),size(mesh.tri)),2 ) ];
  quiver3(centers(:,1),centers(:,2),centers(:,3), ...
        triangledata(:,1),triangledata(:,2),triangledata(:,3),'r');
end

if strcmp( numberelements , 'yes' )
  centers=[ mean( reshape(mesh.xyz(mesh.tri(:),1),size(mesh.tri)),2 ) ...
            mean( reshape(mesh.xyz(mesh.tri(:),2),size(mesh.tri)),2 ) ...
            mean( reshape(mesh.xyz(mesh.tri(:),3),size(mesh.tri)),2 ) ];
  for e=elements
      text( centers(e,1), centers(e,2), centers(e,3) , ...
      sprintf('%d',e),'Color',[0.8 0.5 0.1],...
      'HorizontalAlignment', 'center','VerticalAlignment', 'middle','FontWeight','bold');
  end
end

if strcmp( numberpoints , 'yes' )
  points_with_text= zeros( size(mesh.xyz,1),1);
  for e= elements(:)'
    for ee=1:3
      p= mesh.tri(e,ee);
      if points_with_text(p) == 0
        text( mesh.xyz(p,1) , mesh.xyz(p,2) , mesh.xyz(p,3) , ...
          sprintf('%d',p),'Color',[0 0 0],...
          'HorizontalAlignment', 'left','VerticalAlignment', 'bottom','FontWeight','bold');
        points_with_text(p) = 1;
      end
    end
  end
end

fields= fieldnames( mesh );
for f=1:size(fields,1)
  field= fields{f};
  if ( strncmp( field, 'lmk',3) )
    plot3( mesh.(field)(:,1) , mesh.(field)(:,2) , mesh.(field)(:,3) , lmkls(lmkid,:) , ...
      'MarkerSize',15 ,'displayname', field);
    axis tight;
    lmkid= lmkid+1;
  end
end

if ~bloquea
  if     max( mesh.xyz(:,3) )-min( mesh.xyz(:,3) ) < 1e-6
            view(0,90);
            rotate3d off;
  elseif max( mesh.xyz(:,1) )-min( mesh.xyz(:,1) ) < 1e-6
            view(90,0);
            rotate3d off;
  elseif max( mesh.xyz(:,2) )-min( mesh.xyz(:,2) ) < 1e-6
            view(0,0);
            rotate3d off;
  else
            view(3);
%             rotate3d on;
  end
  hold off;
end

drawnow


% function centers= centertriangles( m )
%   centers(:,1)= mean([ m.xyz(m.tri(:,1),1) m.xyz(m.tri(:,2),1) m.xyz(m.tri(:,3),1) ],2);
%   centers(:,2)= mean([ m.xyz(m.tri(:,1),2) m.xyz(m.tri(:,2),2) m.xyz(m.tri(:,3),2) ],2);
%   centers(:,3)= mean([ m.xyz(m.tri(:,1),3) m.xyz(m.tri(:,2),3) m.xyz(m.tri(:,3),3) ],2);
end


