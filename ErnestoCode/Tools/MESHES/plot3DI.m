function plot3DI( X,Y,Z,D )

if nargin < 2
  D= X;
  [ni,nj,nk] = size(D);
  [X,Y,Z]= ndgrid( 1:ni , 1:nj , 1:nk );
end

for k=1:size(D,3)
  DD(:,:,k)= D(:,:,k)';
  XX(:,:,k)= X(:,:,k)';
  YY(:,:,k)= Y(:,:,k)';
  ZZ(:,:,k)= Z(:,:,k)';
end

mini= XX(1); maxi= XX(end);
minj= YY(1); maxj= YY(end);
mink= ZZ(1); maxk= ZZ(end);

h = slice( XX,YY,ZZ,DD , (mini+maxi)/2, (minj+maxj)/2, (mink+maxk)/2 );

set(h,'FaceColor','interp','EdgeColor','none');

colormap( jet(256) );

lighting flat;

view(3);
axis equal;
xlabel('X');
ylabel('Y');
zlabel('Z');

drawnow;
