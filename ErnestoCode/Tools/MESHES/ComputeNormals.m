function normal = ComputeNormals( mesh )
%
%   normal = ComputeNormals( mesh );
%

e21= mesh.xyz( mesh.tri(:,2),:) - mesh.xyz( mesh.tri(:,1),:);
e31= mesh.xyz( mesh.tri(:,3),:) - mesh.xyz( mesh.tri(:,1),:);

normal = cross( e21 , e31 , 2 );
        
normal = bsxfun(@rdivide,  normal , sqrt( sum( normal.^2 , 2) ) );


%   function M = checkFaces( M )
%     M.tri( any(M.tri == 0,2) | M.tri(:,1) == M.tri(:,2) | M.tri(:,1) == M.tri(:,3) | M.tri(:,2) == M.tri(:,3)  , : ) = [];
%     %flips = any( vtkComputeNormals(M,'SetComputePointNormals',false,'SetComputeCellNormals',true) ./ ComputeNormals( M ) < 0 , 2 );
%     flips = sum( vtkComputeNormals(M,'SetComputePointNormals',false,'SetComputeCellNormals',true) .* ComputeNormals( M ) , 2 ) < 0;
%     M.tri(flips,:) = M.tri(flips,[1 3 2]);
%   end
