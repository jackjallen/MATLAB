function M = CleanMesh( M , threshold )
% 
% mesh= CleanMesh( M )
% 

  if nargin < 2
    threshold = eps(single(1));
  end

  M = vtkCleanPolyData( M , 'SetAbsoluteTolerance', 0 , 'PointMergingOn',[],'ConvertStripsToPolysOn',[],'ConvertPolysToLinesOn',[],'ConvertLinesToPointsOn',[],'ToleranceIsAbsoluteOn',[] );

%   M.tri( any( M.tri == 0 , 2 ) | M.tri(:,1) == M.tri(:,2) | M.tri(:,2) == M.tri(:,3) | M.tri(:,3) == M.tri(:,1) , : ) = [];

%   % Removing zero area triangles
%   e21= M.xyz( M.tri(:,2),:) - M.xyz( M.tri(:,1),:);
%   e31= M.xyz( M.tri(:,3),:) - M.xyz( M.tri(:,1),:);
% %   normal= [ e21(:,2).*e31(:,3)-e21(:,3).*e31(:,2) ...
% %             e21(:,3).*e31(:,1)-e21(:,1).*e31(:,3) ...
% %             e21(:,1).*e31(:,2)-e21(:,2).*e31(:,1) ];
%   normal= e21(:,[2 3 1]).*e31(:,[3 1 2]) - e21(:,[3 1 2]).*e31(:,[2 3 1]);
%   M.tri= M.tri( any( normal , 2 ) , :);

%   E = [ M.tri(:,[1 2]) ; M.tri(:,[2 3]) ; M.tri(:,[3 1]) ];
%   try
%     E = unique([min(E,[],2) max(E,[],2)],'rows');
%   end

%   try,
%     El = sum( ( M.xyz(E(:,1),:) - M.xyz(E(:,2),:) ).^2 , 2 );
%     for id = find( El < threshold^2 )
%       M.tri( ismembc( M.tri , E(id,2) ) ) = E(id,1);
%     end
%   end

%   M.tri( any( M.tri == 0 , 2 ) | M.tri(:,1) == M.tri(:,2) | M.tri(:,2) == M.tri(:,3) | M.tri(:,3) == M.tri(:,1) , : ) = [];

  try
    M = DeletePoints( M , setdiff( 1:size(M.xyz,1) , unique( M.tri(:) ) ) );
  end
  
%   try
%     nM = vtkFeatureEdges( M , 'BoundaryEdgesOff',[],'FeatureEdgesOff',[],'ColoringOff',[],'ManifoldEdgesOff',[],'NonManifoldEdgesOn',[] );
% 
%     if isfield( nM , 'xyz' )  && ~isempty( nM.xyz )
%       [id,p,d]= vtkClosestPoint( M , nM.xyz );
%       id = id( d == 0 );
% 
%       M = DeletePoints( M , id );
%     end
%   end  
  

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Removing auxiliar structures
  if( isfield(M,'PLOC') )
    M= rmfield(M,'PLOC');
  end
  if( isfield(M,'ELOC') )
    M= rmfield(M,'ELOC');
  end
  if( isfield(M,'SCAM') )
    M= rmfield(M,'SCAM');
  end
  if( isfield(M,'ESUP') )
    M= rmfield(M,'ESUP');
  end
  if( isfield(M,'PSUP') )
    M= rmfield(M,'PSUP');
  end









%   for it=1:5
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   % Structured Cartesian Auxiliar Mesh for spatial search and the PointLocator
%   M= CreateLocators( M );
% 
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   % Merging points
%   M= MergePoints( M , threshold );
%   end
%   
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   % Removing auxiliar structures
%   if( isfield(M,'PLOC') )
%     M= rmfield(M,'PLOC');
%   end
%   if( isfield(M,'ELOC') )
%     M= rmfield(M,'ELOC');
%   end
%   if( isfield(M,'SCAM') )
%     M= rmfield(M,'SCAM');
%   end
%   if( isfield(M,'ESUP') )
%     M= rmfield(M,'ESUP');
%   end
%   if( isfield(M,'PSUP') )
%     M= rmfield(M,'PSUP');
%   end
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   % Removing zero area triangles
%   e21= M.xyz( M.tri(:,2),:) - M.xyz( M.tri(:,1),:);
%   e31= M.xyz( M.tri(:,3),:) - M.xyz( M.tri(:,1),:);
%   normal= [ e21(:,2).*e31(:,3)-e21(:,3).*e31(:,2) ...
%             e21(:,3).*e31(:,1)-e21(:,1).*e31(:,3) ...
%             e21(:,1).*e31(:,2)-e21(:,2).*e31(:,1) ];
%   areas= sqrt( sum( normal.^2 , 2));        
%   M.tri= M.tri( areas>0 , :);
% 
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   % Removing unused points
%   n_xyz= size( M.xyz , 1);
%   nonusedpoints= setdiff( 1:n_xyz , unique( M.tri ) );
%   M= DeletePoints( M , nonusedpoints );
end
