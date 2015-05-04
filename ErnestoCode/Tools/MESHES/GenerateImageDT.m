function DT= GenerateImageDT( m , x , y , z , offset )
%
% DT= GenerateImageDT( m , x , y , z)
%

if nargin < 5
  offset= Inf;
end

p= [ x(:) y(:) z(:) ];

bb= BBMesh( m );

ids=  (x(:)>bb(1,1)-offset) & (x(:)<bb(2,1)+offset) ...
    & (y(:)>bb(1,2)-offset) & (y(:)<bb(2,2)+offset) ...
    & (z(:)>bb(1,3)-offset) & (z(:)<bb(2,3)+offset);

DT= x*0; DT(:)= NaN;

[elements, points, dist ]= vtkClosestElement( m , p(ids,:) );

if isfield( m , 'triNORMALS' )
  side= sign( sum( ( points-p(ids,:) ).*m.triNORMALS(elements,:) , 2) );
  if side(1) < 0
    side= -side;
  end
  dist= dist.*side;
end

DT(ids)= dist;
