%%Ernesto's code
%%
fid = fopen('DATA/Labels.csv','r');
SUBJECTS = textscan(fid,'%s %s %s','HeaderLines',1,'Delimiter',',');
fclose(fid);
SUBJECTS = [ SUBJECTS{:} ];

%%

fid = fopen('DATA/TriangleFaces.csv','r');
M.tri = double(cell2mat( textscan(fid,'%d %d %d','Delimiter',',') ));
fclose(fid);

fid = fopen( fullfile('DATA','output-stacom-clean',[SUBJECTS{1,3},'.ED.endo.vertices.csv']) , 'r' );
M.xyz = double(cell2mat( textscan(fid,'%f %f %f','Delimiter',',') ));
fclose(fid);

EL = [ M.xyz(M.tri(:,2),:) - M.xyz(M.tri(:,1),:) ;
       M.xyz(M.tri(:,3),:) - M.xyz(M.tri(:,1),:) ;
       M.xyz(M.tri(:,3),:) - M.xyz(M.tri(:,2),:) ];
EL = sqrt( sum( EL.^2 , 2 ) );
semilogy( sort( EL ) ,'.-' )

% M_clean = vtkCleanPolyData( M ,'SetAbsoluteTolerance',0.1,'SetToleranceIsAbsolute',true,'SetConvertLinesToPoints',false,'SetConvertPolysToLines',false,'SetPointMerging',true,'SetConvertStripsToPolys',true);
% M_clean.tri( any( M_clean.tri == 0 , 2 ) , : ) = [];

M_clean = TidyMesh( M , 0.1 );
M_clean = FixNormals( M_clean );
% plotMESH( M_clean , 'td' , ComputeNormals( M_clean ) )


[NODE_ids,NODE_order] = sort( vtkClosestPoint( M , M_clean.xyz ) );
NODE_iorder( NODE_order ) = 1:numel( NODE_order );
M_clean.tri = NODE_iorder( M_clean.tri );
M_clean.xyz = M.xyz(NODE_ids,:);
% plotMESH( M_clean , 'td' , ComputeNormals( M_clean ) )

B = vtkFeatureEdges( M_clean ,'SetBoundaryEdges',true,'SetFeatureEdges',false,'SetNonManifoldEdges',true,'SetManifoldEdges',false,'SetColoring',false);
B_ids = vtkClosestPoint( M_clean , B.xyz );
% B.tri = B_ids( B.tri(:,1:2) );

M_clean.xyzIsBoundary = -ones( size(M_clean.xyz,1) , 1 );
M_clean.xyzIsBoundary( B_ids ) = 1;

M_clean.triRegion = ones( size(M_clean.tri,1) , 1 );

% NODE_remove = setdiff( 1:size(M.xyz,1) ,  NODE_ids );

%%

for s = 1:size(SUBJECTS,1)
  disp(s)
  
  %{   1 -> ED-ENDO    2 -> ED-EPI    3 -> ES-ENDO   4 -> ES-EPI   %}
  ED_ENDO = readSubject( fullfile('DATA','output-stacom-clean',[SUBJECTS{s,3},'.ED.endo.vertices.csv']) , 1 );
  ED_EPI  = readSubject( fullfile('DATA','output-stacom-clean',[SUBJECTS{s,3},'.ED.epi.vertices.csv'])  , 2 );
  ES_ENDO = readSubject( fullfile('DATA','output-stacom-clean',[SUBJECTS{s,3},'.ES.endo.vertices.csv']) , 3 );
  ES_EPI  = readSubject( fullfile('DATA','output-stacom-clean',[SUBJECTS{s,3},'.ES.epi.vertices.csv'])  , 4 );

  M = AppendMeshes( AppendMeshes( ES_EPI , ES_ENDO ), AppendMeshes( ED_EPI , ED_ENDO ));
  SUBJECTS{s,4} = M;
  
%   plotMESH( closeBloodPool( getPiece( M , 1 ) ) );
%   plotMESH( closeBloodPool( getPiece( M , 2 ) ) );
%   plotMESH( closeBloodPool( getPiece( M , 3 ) ) );
%   plotMESH( closeBloodPool( getPiece( M , 4 ) ) );
%   plotMESH( closeMyocardium( getPiece( M , [1 2] ) ) );
%   plotMESH( closeMyocardium( getPiece( M , [3 4] ) ) );
end

%%

ratio_myocardial_VOL_ES_ED = arrayfun(@(s)MeshVolume( closeMyocardium( getPiece( SUBJECTS{s,4} , [3 4] ) ) ) , vec( 1:size(SUBJECTS,1) ) ) ./...
                             arrayfun(@(s)MeshVolume( closeMyocardium( getPiece( SUBJECTS{s,4} , [1 2] ) ) ) , vec( 1:size(SUBJECTS,1) ) );

EF = arrayfun(@(s)MeshVolume( closeBloodPool( getPiece( SUBJECTS{s,4} , 3 ) ) ) , vec( 1:size(SUBJECTS,1) ) ) ./...
     arrayfun(@(s)MeshVolume( closeBloodPool( getPiece( SUBJECTS{s,4} , 1 ) ) ) , vec( 1:size(SUBJECTS,1) ) );


%%

  arrayfun(@(s)MeshVolume( closeMyocardium( getPiece( SUBJECTS{s,4} , [3 4] ) ) ) , vec( 1:size(SUBJECTS,1) ) )
find(  arrayfun(@(s)MeshVolume( closeMyocardium( getPiece( SUBJECTS{s,4} , [1 2] ) ) ) , vec( 1:size(SUBJECTS,1) ) ) < 0 )


M = getPiece( SUBJECTS{118,4} , [3 4] ); plotMESH( M , 'td' , M.triRegion ); colormap jet; caxis([0 6])
M = getPiece( SUBJECTS{171,4} , [3 4] ); plotMESH( M , 'td' , M.triRegion ); colormap jet; caxis([0 6])
M = getPiece( SUBJECTS{161,4} , [1 2] ); plotMESH( M , 'td' , M.triRegion ); colormap jet; caxis([0 6])
M = getPiece( SUBJECTS{166,4} , [1 2] ); plotMESH( M , 'td' , M.triRegion ); colormap jet; caxis([0 6])

  
   


%% Using Jack Code

fid = fopen('DATA/Labels.csv','r');
SUBJECTS = textscan(fid,'%s %s %s','HeaderLines',1,'Delimiter',',');
fclose(fid);
SUBJECTS = [ SUBJECTS{:} ];

fid = fopen('DATA/TriangleFaces.csv','r');
tri = double(cell2mat( textscan(fid,'%d %d %d','Delimiter',',') ));
fclose(fid);

data = struct();
for s = 1:numel(SUBJECTS)
  fid = fopen( fullfile('DATA','output-stacom-clean',[SUBJECTS{s,3},'.ED.endo.vertices.csv']) , 'r' ); xyz = double(cell2mat( textscan(fid,'%f %f %f','Delimiter',',') )); fclose(fid);
  data(s).diastolic.endo.xyz = xyz;   data(s).diastolic.endo.tri = tri;
  
  fid = fopen( fullfile('DATA','output-stacom-clean',[SUBJECTS{s,3},'.ED.epi.vertices.csv']) , 'r' ); xyz = double(cell2mat( textscan(fid,'%f %f %f','Delimiter',',') )); fclose(fid);
  data(s).diastolic.epi.xyz = xyz;    data(s).diastolic.epi.tri = tri;
  
  fid = fopen( fullfile('DATA','output-stacom-clean',[SUBJECTS{s,3},'.ES.endo.vertices.csv']) , 'r' ); xyz = double(cell2mat( textscan(fid,'%f %f %f','Delimiter',',') )); fclose(fid);
  data(s).systolic.endo.xyz = xyz;    data(s).systolic.endo.tri = tri;
  
  fid = fopen( fullfile('DATA','output-stacom-clean',[SUBJECTS{s,3},'.ES.epi.vertices.csv']) , 'r' ); xyz = double(cell2mat( textscan(fid,'%f %f %f','Delimiter',',') )); fclose(fid);
  data(s).systolic.epi.xyz = xyz;     data(s).systolic.epi.tri = tri;
end





%% Trying to build the bull-eye plots

% EPI = getPiece( SUBJECTS{1,4} , 2 );
% EPI.xyz = reshape(DETERMINEmeanSysEpi,[1089 3]);
dia_myo_mean = reshape(dia_myo_mean, [2178 , 3]);
dia_epi_mean = dia_myo_mean(1090:end,:)
dia_endo_mean = dia_myo_mean(1:1089,:)
hold on
plot3D(dia_endo_mean,'o')
plot3D(dia_epi_mean,'o')
legend

EPI.xyz = dia_epi_mean
EPI.tri = data(1).diastolic.endo.tri;
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
ENDO.xyz = dia_endo_mean;
ENDO.tri = EPI.tri;
ENDO = TidyMesh(ENDO)
% ENDO = data(1).diastolic.endo
[~,~,dia_dEPI2ENDO] = vtkClosestElement( ENDO , EPI.xyz )
% [~,~,dia_dENDO2EPI] = vtkClosestElement( EPI , ENDO.xyz )


figure; subplot 121
hold on
patch('vertices',EPI.xyz, 'faces', EPI.tri,'facecolor','interp','cdata',dia_dEPI2ENDO)
% patch('vertices',ENDO.xyz, 'faces', ENDO.tri,'facecolor','interp','cdata',dia_dENDO2EPI)
view(3)

for m = 1:32
    m
plot3D(EPI.xyz(meridian(m,:),:),'-o','markersize',10)
pause
end

%%
R     = linspace( 0 , 1    , size( meridian , 2 )   );
THETA = linspace( 0 , 2*pi , size( meridian , 1 )+1 );
X     = bsxfun( @times , R , cos(THETA).' );
Y     = bsxfun( @times , R , sin(THETA).' );

 v = dia_dEPI2ENDO( meridian )

% DETERMINEdiaMeanMyoT =  mean(dia_dEPI2ENDOs(data(1).DETERMINE_indices,:),1) ;
% MESAdiaMeanMyoT =  mean(dia_dEPI2ENDOs(data(1).MESA_indices,:),1) ;
% v = DETERMINEdiaMeanMyoT( meridian )
% v = MESAdiaMeanMyoT( meridian )

figure
subplot 122
subplot 221
subplot 222
surf( X , Y , v([1:end 1],:) ,'facecolor','interp'); view(2)
axis equal
cmin = min(min(v([1:end 1],:)))
cmax = max(max(v([1:end 1],:)))
caxis ([ cmin cmax])
title 'DETERMINE mean systole myocardium thickness'
title 'MESA mean systole myocardium thickness'


DETERMINEsysMeanMyoT =  mean(sys_dEPI2ENDOs(data(1).DETERMINE_indices,:),1) ;
MESAsysMeanMyoT =  mean(sys_dEPI2ENDOs(data(1).MESA_indices,:),1) ;

v = DETERMINEsysMeanMyoT( meridian )
v = MESAsysMeanMyoT( meridian )
figure
subplot 223
subplot 224
surf( X , Y , v([1:end 1],:) ,'facecolor','interp'); view(2)
axis equal
cmin = min(min(v([1:end 1],:)))
cmax = max(max(v([1:end 1],:)))
caxis ([ -14 14])
title 'DETERMINE mean systole myocardium thickness'
title 'MESA mean systole myocardium thickness'
