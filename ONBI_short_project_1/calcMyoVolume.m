function[data] = calcMyoVolume(data, myoB)

% for i = 1:400 %all patients
for i = 1:401
    % Find endo and epi boundary points (B.xyz)
    % vtkCleanPolyData(EPI_ED) fix the possible replicated nodes and spurious
    % edges.
    data(i).diastolic.endo.B = vtkFeatureEdges( vtkCleanPolyData(data(i).diastolic.endo) , 'BoundaryEdgesOn' , [] , 'FeatureEdgesOff' , [] );  %extracting the boundary
    data(i).diastolic.endo.B.xyz = data(i).diastolic.endo.B.xyz( [2 1 3:end], : );  %fixing the connectivity.  
    data(i).diastolic.epi.B = vtkFeatureEdges( vtkCleanPolyData(data(i).diastolic.epi) , 'BoundaryEdgesOn' , [] , 'FeatureEdgesOff' , [] );  %extracting the boundary
    data(i).diastolic.epi.B.xyz = data(i).diastolic.epi.B.xyz( [2 1 3:end], : );  %fixing the connectivity.
    
    data(i).systolic.endo.B = vtkFeatureEdges( vtkCleanPolyData(data(i).systolic.endo) , 'BoundaryEdgesOn' , [] , 'FeatureEdgesOff' , [] );  %extracting the boundary
    data(i).systolic.endo.B.xyz = data(i).systolic.endo.B.xyz( [2 1 3:end], : );  %fixing the connectivity.   
    data(i).systolic.epi.B = vtkFeatureEdges( vtkCleanPolyData(data(i).systolic.epi) , 'BoundaryEdgesOn' , [] , 'FeatureEdgesOff' , [] );  %extracting the boundary
    data(i).systolic.epi.B.xyz = data(i).systolic.epi.B.xyz( [2 1 3:end], : );  %fixing the connectivity.
    
%    pat1_endoB = vtkFeatureEdges( vtkCleanPolyData(pat1_endo) , 'BoundaryEdgesOn' , [] , 'FeatureEdgesOff' , [] );  %extracting the boundary
%     pat1_endoB.xyz = pat1_endoB.xyz( [2 1 3:end], : );  %fixing th
%     
%     plot3D(data(i).diastolic.endo.B.xyz)
%     hold on
    
    %Find full shape boundary points (B.xyz)
    % full = full shape without myo lid
    data(i).diastolic.full.xyz = [data(i).diastolic.endo.xyz ; data(i).diastolic.epi.xyz ];
    data(i).diastolic.full.tri = [data(i).diastolic.endo.tri ; data(i).diastolic.epi.tri + size(data(i).diastolic.endo.xyz,1) ];
    data(i).diastolic.full.B = vtkFeatureEdges( vtkCleanPolyData(data(i).diastolic.full) , 'BoundaryEdgesOn' , [] , 'FeatureEdgesOff' , [] );  %extracting the boundary
    data(i).diastolic.full.B.xyz = data(i).diastolic.full.B.xyz( [2 1 3:end], : );  %fixing the connectivity.
  
%     pat1.full.xyz = [pat1_endo.xyz ; pat1_epi.xyz ];
%     pat1.full.tri = [pat1_endo.tri ; pat1_epi.tri + size(pat1_endo.xyz,1) ];
%     pat1.full.B = vtkFeatureEdges( vtkCleanPolyData(pat1.full) , 'BoundaryEdgesOn' , [] , 'FeatureEdgesOff' , [] );  %extracting the boundary
%     pat1.full.B.xyz =  pat1.full.B.xyz( [2 1 3:end], : );  %fixing the connectivity.
%     
    
    %     
%     plot3D(data(i).diastolic.full.xyz)
     
    data(i).systolic.full.xyz = [data(i).systolic.endo.xyz ; data(i).systolic.epi.xyz ];
    data(i).systolic.full.tri = [data(i).systolic.endo.tri ; data(i).systolic.epi.tri + size(data(i).systolic.endo.xyz,1) ];
    data(i).systolic.full.B = vtkFeatureEdges( vtkCleanPolyData(data(i).systolic.full) , 'BoundaryEdgesOn' , [] , 'FeatureEdgesOff' , [] );  %extracting the boundary
    data(i).systolic.full.B.xyz = data(i).systolic.full.B.xyz( [2 1 3:end], : );  %fixing the connectivity.
    
    %join the list of coordinates for endo and epi to be used to make myo lid
    data(i).diastolic.myo.B.xyz = [data(i).diastolic.endo.B.xyz ; data(i).diastolic.epi.B.xyz]; 
    data(i).systolic.myo.B.xyz = [data(i).systolic.endo.B.xyz ; data(i).systolic.epi.B.xyz];
    
 
    
%     pat1.myo.B.xyz = [data(i).diastolic.endo.B.xyz ; data(i).diastolic.epi.B.xyz];
%     
%     
%      plot3D(data(i).diastolic.myo.B.xyz)
%      
    % find nearest points on endo and epi boundaries (identically positioned, but not connected to main shape) and assign them as the
    % boundary of the lid.
    % first arg = a mesh, second arg = a list of point coordinates.
    data(i).diastolic.myo.B.xyz = data(i).diastolic.full.xyz( vtkClosestPoint( data(i).diastolic.full, data(i).diastolic.myo.B.xyz ) , : ); 
    
    data(i).systolic.myo.B.xyz = data(i).systolic.full.xyz( vtkClosestPoint( data(i).systolic.full , data(i).systolic.myo.B.xyz ) , : );


%  pat1.myo.B.xyz = pat1.full.xyz( vtkClosestPoint(pat1.full, pat1.myo.B.xyz ) , : );
%     

%     plot3D(data(i).diastolic.myo.B.xyz)
    
    %Make fully closed myocardium volume by appending epi surface, endo
    %surface and myo lid.
    % myoB is the vertices for covering the top of the myocardium.
    data(i).diastolic.myo.xyz = [ data(i).diastolic.endo.xyz ; data(i).diastolic.myo.B.xyz ; data(i).diastolic.epi.xyz ];
    data(i).diastolic.myo.tri = [ data(i).diastolic.endo.tri; myoB.tri + size( data(1).diastolic.endo.xyz , 1 ) ; data(i).diastolic.epi.tri + size( data(1).diastolic.endo.xyz , 1 ) + size(data(i).diastolic.myo.B.xyz,1) ];
    
    data(i).systolic.myo.xyz = [ data(i).systolic.endo.xyz ; data(i).systolic.myo.B.xyz ; data(i).systolic.epi.xyz ];
    data(i).systolic.myo.tri = [ data(i).systolic.endo.tri; myoB.tri + size( data(1).systolic.endo.xyz , 1 ) ; data(i).systolic.epi.tri + size( data(1).systolic.endo.xyz , 1 ) + size(data(i).systolic.myo.B.xyz,1) ];
     

%         pat1.myo.xyz = [ pat1_endo.xyz ; pat1.myo.B.xyz ; pat1_epi.xyz ];
% 
%      pat1.myo.tri = [ pat1_endo.tri; myoB.tri + size( pat1_endo.xyz , 1 ) ; pat1_epi.tri + size( pat1_endo.xyz , 1 ) + size(pat1.myo.B.xyz,1) ];
%     
%     plot3D(data(i).diastolic.myo.xyz)
%     hold on
%      plot3D(data(i).systolic.myo.xyz)
%     
%   patch('vertices',data(i).diastolic.myo.xyz,'faces',data(i).diastolic.myo.tri,'FaceColor','r')
%   
    %make sure that the normal to each triangle points outwards.
    data(i).diastolic.myo = FixNormals(data(i).diastolic.myo );
    
    data(i).systolic.myo = FixNormals(data(i).systolic.myo );
    
%      pat1.myo = FixNormals(pat1.myo );
%      
    %calculate volume of myocardium.  
   [data(i).diastolic.myoVolume, data(i).diastolic.myoCenterOfMass] = MeshVolume( data(i).diastolic.myo );    
   [data(i).systolic.myoVolume, data(i).systolic.myoCenterOfMass] = MeshVolume( data(i).systolic.myo );
   
%      [pat1.myoVolume, pat1.myoCenterOfMass] = MeshVolume( pat1.myo );
    
    % %Compare myo volume with epi volume.
    % %transformed_data(i).diastolic.myodifference_volume = prod(diff( BBMesh( transformed_data(i).diastolic.myoB_full ) , 1  , 1 ) ) - transformed_data(i).diastolic.myoVolume ;   %%it shoud be positive!!
    % transformed_data(i).diastolic.myodifference_volume =  transformed_data(i).diastolic.epi.volume - transformed_data(i).diastolic.myoVolume ;   %%it shoud be positive!!
    % transformed_data(i).systolic.myodifference_volume = transformed_data(i).systolic.epi.volume - transformed_data(i).systolic.myoVolume ;   %%it shoud be positive!!
end

end
