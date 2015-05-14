%Read files
vertices = data(1).diastolic.endo.xyz;
trifac = data(1).diastolic.endo.tri;


%Select only the non-duplicated vertices
[vertices,IA,IC] = unique(vertices,'rows');


%Clean the face matrix
trifac=IC(trifac);  %Change the vertex numbers for the new ones
non_degen=find(trifac(:,1)~=trifac(:,2) & trifac(:,1)~=trifac(:,3) & trifac(:,2)~=trifac(:,3)); %Find all non-degenerate triangles, i.e. those in which vertices are not repeated
trifac=trifac(non_degen,:); %Keep only non-degenerate triangles in trifac


%For visualisation only
p.vertices=vertices;
p.faces=trifac;
figure
patch(p,'EdgeColor','black','FaceColor','red');
axis equal; axis tight;