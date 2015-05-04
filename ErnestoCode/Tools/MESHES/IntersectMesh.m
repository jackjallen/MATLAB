function [m, mm1, mm2]= IntersectMesh( m1, m2)
% 
% m= IntersectMesh( m1, m2)
% 

  mm2= m2;
  mm1= m1;
  for i=1:size(m1.tri,1)
    mm2= cut_( mm2, m1.xyz( m1.tri(i,1),:),m1.xyz( m1.tri(i,2),:),m1.xyz( m1.tri(i,3),:) );
  end

  for i=1:size(m2.tri,1)
    mm1= cut_( mm1, m2.xyz( m2.tri(i,1),:),m2.xyz( m2.tri(i,2),:),m2.xyz( m2.tri(i,3),:) );
  end

%   mm1= vtkcut_( m1, m2);
%   mm2= vtkcut_( m2, m1);

%   mm1= CleanMesh(mm1);
%   mm1.triNORMALS= vtkComputeNormals( mm1 );
% 
%   mm2= CleanMesh(mm2);
%   mm2.triNORMALS= vtkComputeNormals( mm2 );

  mmm1= DeletePoints( mm1 , find( vtkIsInside(m2, mm1.xyz) < 0) );
  mmm2= DeletePoints( mm2 , find( vtkIsInside(m1, mm2.xyz) < 0) );
  
  m.xyz= [ mmm1.xyz ; mmm2.xyz ];
  m.tri= [ mmm1.tri ; mmm2.tri + size(mmm1.xyz,1) ];
  m.triNORMALS= [ mmm1.triNORMALS ; mmm2.triNORMALS ];


function m= cut_( m, p1, p2, p3)
  bbs= [ m.xyz(m.tri(:,1),:) m.xyz(m.tri(:,2),:) m.xyz(m.tri(:,3),:) ];
  bbs= [ min(bbs(:,[1 4 7]),[],2)  max(bbs(:,[1 4 7]),[],2) ...
         min(bbs(:,[2 5 8]),[],2)  max(bbs(:,[2 5 8]),[],2) ...
         min(bbs(:,[3 6 9]),[],2)  max(bbs(:,[3 6 9]),[],2) ];

  bb = [ min([p1(1) p2(1) p3(1)]) max([p1(1) p2(1) p3(1)]) ...
         min([p1(2) p2(2) p3(2)]) max([p1(2) p2(2) p3(2)]) ...
         min([p1(3) p2(3) p3(3)]) max([p1(3) p2(3) p3(3)]) ];
       
  p12= p2-p1;
  p13= p3-p1;

  zz = cross( p12, p13);
  yy = cross( zz , p12);
  xx = p12;

  xx= xx/norm(xx)+p1;
  yy= yy/norm(yy)+p1;
  zz= zz/norm(zz)+p1;

  OF= [ p1 1; xx 1; yy 1; zz 1]';
  NF= [ 0 0 0 1; 1 0 0 1; 0 1 0 1; 0 0 1 1]';
  
  T= NF*inv(OF);
  
  xyz= m.xyz'; xyz(4,:)=1;
  xyz= T*xyz ; xyz= xyz(1:3,:)./repmat(xyz(4,:),3,1);
  m.xyz= xyz';
  
  z= [ m.xyz(m.tri(:,1),3) m.xyz(m.tri(:,2),3) m.xyz(m.tri(:,3),3) ];
  t_to_cut= find( abs( sum(sign(z),2) ) < 3 );
  
  for t=t_to_cut'
    bbi= [ max([bb(1) bbs(t,1)])  min([bb(2) bbs(t,2)]) ...
           max([bb(3) bbs(t,3)])  min([bb(4) bbs(t,4)]) ...
           max([bb(5) bbs(t,5)])  min([bb(6) bbs(t,6)]) ];
         
    if all( (bbi(2:2:end)-bbi(1:2:end)) > 0 )

      if sign( z(t,1) )== sign( z(t,2) )
        f= 0.5;
      else
        f = abs( z(t,2)/(z(t,1)-z(t,2)) );
      end
      m12= f*m.xyz(m.tri(t,1),:) + (1-f)*m.xyz(m.tri(t,2),:);

      if sign( z(t,2) )== sign( z(t,3) )
        f= 0.5;
      else
        f = abs( z(t,3)/(z(t,2)-z(t,3)) );
      end
      m23= f*m.xyz(m.tri(t,2),:) + (1-f)*m.xyz(m.tri(t,3),:);

      if sign( z(t,3) )== sign( z(t,1) )
        f= 0.5;
      else
        f = abs( z(t,1)/(z(t,3)-z(t,1)) );
      end
      m31= f*m.xyz(m.tri(t,3),:) + (1-f)*m.xyz(m.tri(t,1),:);

      np= size( m.xyz, 1);

      m.xyz= [ m.xyz ; m12 ; m23 ; m31 ];
      m.tri= [m.tri ; [m.tri(t,2) np+1 np+2]; [m.tri(t,3) np+2 np+3]; [np+1 np+2 np+3] ];
      m.triNORMALS= [m.triNORMALS ; repmat( m.triNORMALS(t,:) ,3,1) ];
      m.tri(t,2:3)= [ np+1 np+3];
    end
  end

  xyz= m.xyz'; xyz(4,:)=1;
  xyz= inv(T)*xyz ; xyz= xyz(1:3,:)./repmat(xyz(4,:),3,1);
  m.xyz= xyz';
