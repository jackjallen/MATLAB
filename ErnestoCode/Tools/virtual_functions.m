function varargout = virtual_functions( fn , varargin )
[ varargout{1:nargout} ] = feval( fn , varargin{:} ); 
end

%VFUNCTION  IPfun                           %created: 11-Nov-2011 17:06:41
function IPfun(ip)
  %variables by vfunction
  h_Y                            = 181.001953125;
  h_X                            = 180.001953125;
  h_T                            = 185.001953125;
  h_STREAM_ends                  = 244.001953125;
  h_STREAM                       = 243.001953125;
  h_QUIVER                       = 239.001953125;
  h_P                            = 184.001953125;
  h_N                            = 182.001953125;
  Y                              = [-1;-0.1;-2.9;-2.4;2.2];
  M                              = [0.3,1;-0.4,1;-2,1;-0.9,1;-1.1,1];
  %END_variables by vfunction
 
set( [ h_STREAM h_STREAM_ends ] , 'visible','off' );
set( h_QUIVER , 'visible','on'  );
ab = squeeze( linspacen( ip.xyz(1,1:2) , ip.xyz(2,1:2) , 100 ) ).';
dirs = ab*0;
for i = 1:size(ab,1)
dirs(i,:) = vect( LAD_dir( M , Y , vec( ab(i,:) ) ) );
end
set( h_QUIVER , 'xdata' , ab(:,1) , 'ydata' , ab(:,2) , 'udata' , dirs(:,1) , 'vdata' , dirs(:,2) );
set( h_X      , 'ydata' , dirs(:,1) );
set( h_Y      , 'ydata' , dirs(:,2) );
D = normalize( ip.xyz(2,1:2) - ip.xyz(1,1:2) );
T = dirs*[ D.' , [D(2);-D(1)] ];
set( h_T      , 'ydata' , T(:,1) );
set( h_P      , 'ydata' , T(:,2) );
set( h_N      , 'ydata' , sqrt(sum(dirs.^2,2)));
end


%VFUNCTION  IPstream                        %created: 11-Nov-2011 17:06:50
function IPstream(ip)
  %variables by vfunction
  h_STREAM_ends                  = 244.001953125;
  h_STREAM                       = 243.001953125;
  h_QUIVER                       = 239.001953125;
  Y                              = [-1;-0.1;-2.9;-2.4;2.2];
  M                              = [0.3,1;-0.4,1;-2,1;-0.9,1;-1.1,1];
  %END_variables by vfunction
 
set( [ h_STREAM h_STREAM_ends ], 'visible','on'  );
set( h_QUIVER , 'visible','off' );
ab = squeeze( linspacen( ip.xyz(1,1:2) , ip.xyz(2,1:2) , 20 ) ).';
st = []; st_end = [];
for i = 1:size(ab,1)
[T,this_st] = ode45( @(t,x) LAD_dir(M,Y,x) ,[0 2],ab(i,:) ,...
odeset('RelTol',1e-10,'AbsTol',1e-10,'InitialStep',1e-8) );
st     = [ st ; this_st ; NaN NaN ];
st_end = [ st_end ; this_st(end,:) ];
set( h_STREAM , 'xdata', st(:,1) , 'ydata',st(:,2) );
set( h_STREAM_ends , 'xdata', st_end(:,1) , 'ydata',st_end(:,2) );
drawnow;
end
end


%VFUNCTION  VF1                             %created: 14-Nov-2011 10:51:00
function [a,b]=VF1(c,d)
  %variables by vfunction
  r                              = 6;
  a                              = 123456789;
  %END_variables by vfunction
 
aa = 7;   a = a + c*d+r; r = 4;
b = d+c;
b = b*5 + aa;
end


%VFUNCTION  rigid                           %created: 27-Feb-2012 15:29:34
function dy=rigid(t,y)
  %variables by vfunction
  %END_variables by vfunction
 
dy = zeros(3,1);
dy(1) = y(2) * y(3);
dy(2) = -y(1) * y(3);
dy(3) = -0.51 * y(1) * y(2);
end


%VFUNCTION  setPoint                        %created: 08-Mar-2012 14:46:41
function setPoint( hP , hp , M , S , hp2 )
  %variables by vfunction
  %END_variables by vfunction
 
if isempty( pressedkeys )
xyz = IntersectSurfaceRay( hP );
set( hp , 'XData' , xyz(1) , 'YData', xyz(2) , 'ZData', xyz(3) );
[tri,kk,kk,bary] = vtkClosestElement( M , xyz );
xyz = Barycentric2XYZ( S , tri , bary );
xyz = normalize( xyz )*100;
set( hp2 , 'XData' , xyz(1) , 'YData', xyz(2) , 'ZData', xyz(3) );
end
end


%VFUNCTION  settingPoint                    %created: 09-Mar-2012 15:47:16
function settingPoint( R1 , P1 , R2 , P2 )
  %variables by vfunction
  hpSULC                         = 179.018310546875;
  hpCURV                         = 188.018310546875;
  hp2                            = 196.018310546875;
  hp1                            = 192.018310546875;
  hS                             = 176.018310546875;
  %END_variables by vfunction
 
if ~isempty( pressedkeys ), return; end
xyz = IntersectSurfaceRay( hS );
set( [hpSULC hpCURV] , 'XData' , xyz(1) , 'YData', xyz(2) , 'ZData', xyz(3) );
[tri,kk,kk,bary] = vtkClosestElement( R1 , xyz );
nxyz = Barycentric2XYZ( P1 , tri , bary );
set( hp1 , 'XData' , nxyz(1) , 'YData', nxyz(2) , 'ZData', nxyz(3) );
[tri,kk,kk,bary] = vtkClosestElement( R2 , xyz );
nxyz = Barycentric2XYZ( P2 , tri , bary );
set( hp2 , 'XData' , nxyz(1) , 'YData', nxyz(2) , 'ZData', nxyz(3) );
end


%VFUNCTION  matchPoints                     %created: 12-Mar-2012 23:07:21
function [E,D]=matchPoints( V , X , Y , LLfft , Kfft , lambda )
  %variables by vfunction
  kk                             = evalin('base','kk');
  %END_variables by vfunction
 
X1 = squeeze( expm( V , X ) );
Em = sum( ( X1(:) - Y(:) ).^2 );
LLV = real( ifftn( bsxfun(@times, fftn(V.data) , LLfft ) ) );
Er = LLV(:).' * V.data(:);
vv = voxelvolume( V );
E = vv * Er + lambda * Em;
if nargout > 1
FORCES = 2*( X1 - Y );
[kk,D] = dEvolvePointsOn3DGrid( V.data , V.X , V.Y , '2d' , X ,'forces', FORCES,'circular','maxSteps',100,'jac_times' );
D = lambda * D + 2*vv*LLV;
if ~isequal( Kfft , 1 )
D = real( ifftn( bsxfun( @times, fftn( D ) , Kfft ) ) );
end
end
end


%VFUNCTION  energy                          %created: 19-Apr-2012 23:51:19
function [E,dE]=energy( X0 )
  %variables by vfunction
  Xs                             = evalin('base','Xs');
  %END_variables by vfunction
 
[E,dE] = Log_ST( X0 , Xs , 'fast' );
E = sum( E.^2 );
dE = sum( dE , 3 );
end


%VFUNCTION  disp_and_blink                  %created: 10-Jun-2012 23:59:34
function disp_and_blink(h),
  %variables by vfunction
  popupmenu                      = 171.0057373046875;
  %END_variables by vfunction
 
set(popupmenu,'value',find(strcmp(get(popupmenu,'string'),get(h,'tag')))); reds = findall(gcf,'Type','Patch','FaceColor',[1 0 0]);  for r= reds(:)',set(r,'Facecolor',get(r,'UserData'));  end;  set(h,'UserData',get(h,'FaceColor'),'FaceColor',[1 0 0]); set(findall(gcf,'Type','Patch','hittest','on'),'Visible','off'); set(h,'visible','on'); pause(.2);set(findall(gcf,'Type','Patch','hittest','on'),'Visible','on');
end


%VFUNCTION  loadBVI                         %created: 25-Jun-2012 10:20:17
function D=loadBVI( fn , mat )
  %variables by vfunction
  %END_variables by vfunction
 
D = loadI3D( fn );
DT = read_bvMINF([fn,'.minf'],false); DT.voxel_size(end:4) = 1;
DT = diag( DT.voxel_size(1:4) ) * DT.storage_to_memory * diag( 1./[ D.deltaX(1) , D.deltaY(1) , D.deltaZ(1) , 1 ] );
D.SpatialTransform = DT;
D = rot90( D , rot90( D , mat ) );
D = D.data;
end


%VFUNCTION  QP_dot                          %created: 10-Jul-2012 10:22:50
function QP_d=QP_dot( QP )
  %variables by vfunction
  S                              = diag(repmat(eye(1,2),[1,2]));
  %END_variables by vfunction
 
Q = tosquare( QP(1:end/2) );
P = tosquare( QP(end/2+1:end) );
U = tosquare( S*vec( Q.'*P ) );
QP_d = [ vec( Q*U ) ; - vec( P*U.' ) ];
end


%VFUNCTION  sum2                            %created: 27-Aug-2012 19:16:43
function x=sum2(x)
  %variables by vfunction
  %END_variables by vfunction
 
while numel(x)>1
if mod(numel(x),2), x(numel(x)+1) = 0; end
x = sum(reshape(x,2,[]),1);
end
end


%VFUNCTION  PLOTEA                          %created: 16-Oct-2012 23:59:42
function PLOTEA(g,t)
  %variables by vfunction
  hg_Y                           = [894.010986328125;895.010498046875;896.010498046875;897.010498046875;898.010498046875;899.010498046875;900.010498046875;901.010498046875;902.010498046875;903.010498046875];
  hg_MODEL                       = 904.010498046875;
  g_to_data                      = evalin('base','g_to_data');
  T                              = evalin('base','T');
  %END_variables by vfunction
 
TS = bsxfun(@plus,T,t);
for tt=1:size(TS,1)
set( hg_Y(tt),'XData',TS(tt,:) );
end
set( hg_MODEL , 'YData' , g_to_data(g) );
end


%VFUNCTION  LIK                             %created: 17-Jul-2013 11:09:56
function lik=LIK( s , l , cf )
  %variables by vfunction
  Y                              = [8.141814707233607;1.758891988150772;5.32781061393832;5.5043018279917035;5.810535602971414;7.620874779016235;6.848916589524103;0.48114525489420895;6.942236999347008;5.321900958348639;4.87377570492371;7.202715824054042;7.257859992067113;8.926170629002819;3.464707037951571];
  Mr                             = evalin('base','Mr');
  Mf                             = [1,0.9501292851471754;1,0.23113851357428783;1,0.6068425835417866;1,0.48598246870929973;1,0.8912989661489016;1,0.7620968330273947;1,0.45646766516834136;1,0.0185036432482244;1,0.8214071642952533;1,0.4447033643531942;1,0.6154323481000947;1,0.7919370374270354;1,0.9218129707448025;1,0.7382072458106653;1,0.176266144494618];
  LtL                            = diag([2:7]);
  %END_variables by vfunction
 
v = zeros(size(Y));
for i = 1:numel(Y),  v(i) = l^2 * Mr(i,:)* LtL \ Mr(i,:).' + s^2; end
lik = sum( log(v) ) + sum( ( ( Mf*cf - Y ).^2 )./v );
end


%VFUNCTION  glyphGL                         %created: 29-Oct-2013 16:31:06
function h=glyphGL( J ,varargin )
  %variables by vfunction
  %END_variables by vfunction
 
g = [0 1 1 -1 -1 1 1;0 0 1 1 -1 -1 0];
vertices = [];
for n = 1:size(J,3)
vertices = [ vertices ; (J(:,:,n)*g).' ];
end; vertices(:,3)=0;
h=patch('vertices',vertices,'faces',reshape(1:size(vertices,1),size(g,2),[]).','facecolor','none',varargin{:});
if ~ishold(get(h,'Parent')), axis( get(h,'Parent') , 'equal' ); end
end


%VFUNCTION  glyphJtJ                        %created: 29-Oct-2013 17:21:02
function h=glyphJtJ( JtJ ,varargin )
  %variables by vfunction
  h                              = 260.0028076171875;
  %END_variables by vfunction
 
g = linspace(0,2*pi,151); g=[cos(g);sin(g)];
vertices = [];
for n = 1:size(JtJ,3)
[R,S]=eig( JtJ(:,:,n) );
vertices = [ vertices ; ( R*diag(1./sqrt(diag(S)))*g ).' ];
end; vertices(:,3)=0;
h=patch('vertices',vertices,'faces',reshape(1:size(vertices,1),size(g,2),[]).','facecolor','none',varargin{:});
if ~ishold(get(h,'Parent')), axis( get(h,'Parent') , 'equal' ); end
end


%VFUNCTION  close_contours                  %created: 29-Mar-2014 19:19:28
function x = close_contours(x),
  %variables by vfunction
  %END_variables by vfunction
 
x(end+1,:) = NaN; x = cell2mat( cellfun( @(c) [ c(1:end-1,:) ; c(1,:) ; c(end,:) ] , mat2cell(x, diff( [ 0 ; find( any(isnan(x),2) ) ] ), 3 ) ,'UniformOutput',false ) );
end


%VFUNCTION  MatchingEnergy                  %created: 23-Jun-2014 14:32:48
function [E,D]=MatchingEnergy( V , X , Y , LLfft , Kfft , lambda )
  %variables by vfunction
  kk                             = evalin('base','kk');
  D                              = evalin('base','D');
  %END_variables by vfunction
 
X1 = squeeze( expm( V , X ) );
Em = sum( ( X1(:) - Y(:) ).^2 );
LLV = real( ifftn( bsxfun(@times, fftn(V.data) , LLfft ) ) );
Er = LLV(:).' * V.data(:);
vv = voxelvolume( V );
E = vv * Er + lambda * Em;
if nargout > 1
FORCES = 2*( X1 - Y );
[kk,D] = dEvolvePointsOn3DGrid( V.data , V.X , V.Y , '2d' , X ,'forces', FORCES,'circular','maxSteps',100,'jac_times' );
D = lambda * D + 2*vv*LLV;
if ~isequal( Kfft , 1 )
D = real( ifftn( bsxfun( @times, fftn( D ) , Kfft ) ) );
end
end
end


%VFUNCTION  drawV                           %created: 23-Jun-2014 14:32:49
function out=drawV( V , h )
  %variables by vfunction
  p                              = 3;
  Y                              = [30,11;25,15;20,5];
  X                              = [10*ones(1,2);10,20;38,24];
  %END_variables by vfunction
 
out = true;
if nargin > 1, figure( h ); clf; end
Qs = expm( V , X , 'flow',[0 0:0.01:1] );
quiver( V )
for p = 1:size(Qs,1)
hplot( X(p,1) , X(p,2) , 'om' ,'markerFaceColor','k','markersize',8);
hplot( squeeze( Qs(p,:,1) ) , squeeze( Qs(p,:,2) ) , '.-r' );
hplot( Y(p,1) , Y(p,2) , 'om' ,'markerFaceColor','m','markersize',8);
end
end


%VFUNCTION  plotsubcortical                 %created: 02-Sep-2014 19:13:43
function hg=plotsubcortical( S )
  %variables by vfunction
  %END_variables by vfunction
 
hg = zeros(14,1);
hg(1)=plotMESH( S.LAccu ,'edgecolor','none','facecolor',[1 0 0]);
hg(2)=plotMESH( S.LAmyg ,'edgecolor','none','facecolor',[0 1 0]);
hg(3)=plotMESH( S.LCaud ,'edgecolor','none','facecolor',[0 0 1]);
hg(4)=plotMESH( S.LHipp ,'edgecolor','none','facecolor',[0.8 0.8 0]);
hg(5)=plotMESH( S.LPall ,'edgecolor','none','facecolor',[1 0 1]);
hg(6)=plotMESH( S.LPuta ,'edgecolor','none','facecolor',[0 .8 .8]);
hg(7)=plotMESH( S.LThal ,'edgecolor','none','facecolor',[0.98 0.6 0.5]);
hg(8)=plotMESH( S.RAccu ,'edgecolor','none','facecolor',[1 0 0]);
hg(9)=plotMESH( S.RAmyg ,'edgecolor','none','facecolor',[0 1 0]);
hg(10)=plotMESH( S.RCaud ,'edgecolor','none','facecolor',[0 0 1]);
hg(11)=plotMESH( S.RHipp ,'edgecolor','none','facecolor',[0.8 0.8 0]);
hg(12)=plotMESH( S.RPall ,'edgecolor','none','facecolor',[1 0 1]);
hg(13)=plotMESH( S.RPuta ,'edgecolor','none','facecolor',[0 .8 .8]);
hg(14)=plotMESH( S.RThal ,'edgecolor','none','facecolor',[0.98 0.6 0.5]);
end


