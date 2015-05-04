function areas= MeshAreaTriangles( M , e )
% 
% a= AreaTriangle( mesh , e )
% 


  if nargin < 2
    e = 1:size(M.tri);
  end
  e = e(:);

  l1 = M.xyz(M.tri(e,2),:) - M.xyz(M.tri(e,1),:);
  l2 = M.xyz(M.tri(e,3),:) - M.xyz(M.tri(e,1),:);

  areas = [ l1(:,2).*l2(:,3) - l1(:,3).*l2(:,2) ...
            l1(:,3).*l2(:,1) - l1(:,1).*l2(:,3) ...
            l1(:,1).*l2(:,2) - l1(:,2).*l2(:,1) ];

  areas = sqrt( sum( areas.^2 , 2 ) ) / 2;

end
