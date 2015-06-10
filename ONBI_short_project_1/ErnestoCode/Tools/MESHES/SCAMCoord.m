function ijk= SCAMCoord( p , mesh )
% 
% ijk= SCAMCoord( points , mesh )
% 

% O= repmat([ mesh.SCAM.Ox mesh.SCAM.Oy mesh.SCAM.Oz ],size(p,1),1); 
% D= repmat([ mesh.SCAM.Dx mesh.SCAM.Dy mesh.SCAM.Dz ],size(p,1),1); 
% 
% ijk= ceil( (p - O)./D );
% ijk( ijk < 1 )= 1;
% ijk( find(ijk(:,1) > mesh.SCAM.Nx),1 )= mesh.SCAM.Nx;
% ijk( find(ijk(:,2) > mesh.SCAM.Ny),2 )= mesh.SCAM.Ny;
% ijk( find(ijk(:,3) > mesh.SCAM.Nz),3 )= mesh.SCAM.Nz;
