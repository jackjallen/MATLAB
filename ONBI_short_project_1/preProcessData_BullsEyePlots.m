
%Ernesto's code

shape = reshape(dia_sys_myo_new_shape_minues, [2*2178 , 3]);
shape = shape(1:2178,:);

dia_epi = shape(1090:end,:) ; 
dia_endo = shape(1:1089,:) ;
hold on
plot3D(dia_endo,'o')
plot3D(dia_epi,'o')
legend

EPI.xyz = dia_epi
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
%%

%%
ENDO.xyz = dia_endo;
ENDO.tri = EPI.tri;
ENDO = TidyMesh(ENDO)
% ENDO = data(1).diastolic.endo
[~,~,dEPI2ENDO] = vtkClosestElement( ENDO , EPI.xyz )
% [~,~,dia_dENDO2EPI] = vtkClosestElement( EPI , ENDO.xyz )

figure; subplot 121
hold on
patch('vertices',EPI.xyz, 'faces', EPI.tri,'facecolor','interp','cdata',dEPI2ENDO)
% patch('vertices',ENDO.xyz, 'faces', ENDO.tri,'facecolor','interp','cdata',dia_dENDO2EPI)
view(3)
% 
% for m = 1:32
%     m
% plot3D(EPI.xyz(meridian(m,:),:),'-o','markersize',10)
% pause
% end

%%
R     = linspace( 0 , 1    , size( meridian , 2 )   );
THETA = linspace( 0 , 2*pi , size( meridian , 1 )+1 );
X     = bsxfun( @times , R , cos(THETA).' );
Y     = bsxfun( @times , R , sin(THETA).' );

v = dEPI2ENDO(meridian);
normaliseddEPI2ENDO = dEPI2ENDO/max(max(v([1:end 1],:)));
v = normaliseddEPI2ENDO( meridian );

figure
subplot 121
patch('vertices',EPI.xyz, 'faces', EPI.tri,'facecolor','interp','cdata',normaliseddEPI2ENDO)
view(3)
axis equal
cmin = min(min(normaliseddEPI2ENDO([1:end 1],:)));
cmax = max(max(normaliseddEPI2ENDO([1:end 1],:)));
caxis ([ cmin cmax])
title 'mean diastole shape for all training shapes'
xlabel 'x'
ylabel 'y'
zlabel 'z'
subplot 122
surf( X , Y , v([1:end 1],:) ,'facecolor','interp'); view(2)
axis ([-1 1 -1 1])
cmin = min(min(v([1:end 1],:)));
cmax = max(max(v([1:end 1],:)));
caxis ([cmin cmax]);
title 'myocardium thickness'
c = colorbar;
c.Label.String = 'normalised thickness';