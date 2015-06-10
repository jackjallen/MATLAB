%% make Meridians
%Ernesto's code

epi = shape(1090:end,:) ; 
endo = shape(1:1089,:) ;
% hold on
% plot3D(endo,'o')
% plot3D(epi,'o')
% legend

EPI.xyz = epi
EPI.tri = data(1).diastolic.endo.tri(1:2112,:)
EPI = TidyMesh(EPI)
valence = accumarray( EPI.tri(:) , 1 );
[~,APEXid] = max(valence);

rings = {APEXid};
while 1
  rings{end+1,1} = setdiff( EPI.tri( any( ismember( EPI.tri , cell2mat(rings) ) , 2 ) , : ) , cell2mat(rings) );
  if isempty( rings{end} ), break; end
end
rings( end , : ) = [];

%reorder the second ring
r2 = min( rings{2} );
rings{2} = setdiff( rings{2} , r2 );
while ~isempty( rings{2} )
    r2(end+1,1) = min( intersect( setdiff( EPI.tri( any( ismember( EPI.tri , r2(end) ) , 2 ) , : ),r2 ) , rings{2} ) );
    rings{2} = setdiff( rings{2} , r2 );
end
rings{2} = r2;

clear meridian
for m = 1:numel( rings{2} )
    
  meridian(m,1) = APEXid;
  meridian(m,2) = rings{2}(m);

   for r = 2:size( rings , 1 )-1
    nodes_in_next_ring = rings{r+1};
    
    distance = sqrt( sum( bsxfun( @minus , EPI.xyz( nodes_in_next_ring , : ) , EPI.xyz( meridian(m,r) , : ) ).^2 , 2 ) );

    %      point2nextRing= sqrt( sum( bsxfun( @minus , EPI.xyz( nodes_in_next_ring , : ) , EPI.xyz( meridian(m,r) , : ) ).^2 , 2 ) );

%      nextpoint2nextRing = sqrt( sum( bsxfun( @minus , EPI.xyz( nodes_in_next_ring , : ) , EPI.xyz( rings{r}(m+1) , : ) ).^2 , 2 ) );
%      distance3 = sqrt( sum( bsxfun( @minus , EPI.xyz( nodes_in_next_ring , : ) , EPI.xyz( meridian(m-1,r) , : ) ).^2 , 2 ) );
% distance = point2nextRing - nextpoint2nextRing;
% distance = abs(distance);

[~,closestID] = min(distance);
    meridian(m,r+1) = (nodes_in_next_ring (closestID ));
     rings{r+1} = setdiff( rings{r+1} , meridian(m,r+1) );

   end
  
 
end