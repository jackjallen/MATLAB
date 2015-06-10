


diaPts = [1 : 1089];
sysPts = [1090 : 2178];

R     = linspace( 0 , 1    , size( meridian , 2 )   );
THETA = linspace( 0 , 2*pi , size( meridian , 1 )+1 );
X     = bsxfun( @times , R , cos(THETA).' );
Y     = bsxfun( @times , R , sin(THETA).' );

%%
factor1=nStd
dia_sys_dEPI2ENDOs_plus(:,visualisingMode1) = dia_sys_dEPI2ENDOs_mean(:) + principle_dia_sys_dEPI2ENDOs_eigenvectors(:,visualisingMode1)*(factor1)*dia_sys_dEPI2ENDOs_max_b(visualisingMode1);
factor2=-nStd
dia_sys_dEPI2ENDOs_minus(:,visualisingMode1) = dia_sys_dEPI2ENDOs_mean(:) + principle_dia_sys_dEPI2ENDOs_eigenvectors(:,visualisingMode1)*(factor2)*dia_sys_dEPI2ENDOs_max_b(visualisingMode1);
visualisingMode1

clear shape
clear shape1
clear shape2
clear shape3


shape1 = dia_sys_dEPI2ENDOs_minus(:,visualisingMode1)';
figure('name', 'PCA on myocardium thicknesses')
 subplot 221
colormap jet
shape = shape1(sysPts)
% run('make_meridians.m');

%%
v = shape(meridian);
surf( X , Y , v([1:end 1],:) ,'facecolor','interp'); view(2)
axis ([-1 1 -1 1])
% cmin = min(min(v([1:end 1],:)));
% cmax = max(max(v([1:end 1],:)));
% caxis ([cmin cmax]);
title (['sys, mode:', num2str(visualisingMode1),' ',num2str(factor2),'std']);
max_v(1) = max(max(v));
min_v(1) = min(min(v));
 caxis ([cmin  cmax ]);
 c = colorbar;
%  c.Label.String = 'thickness (mm)';
  axis square
 
 shape2 =dia_sys_dEPI2ENDOs_mean';
% figure('name', 'PCA on myocardium thicknesses')
 subplot 222
colormap jet
shape = shape2(sysPts)
% run('make_meridians.m');
v = shape(meridian);
surf( X , Y , v([1:end 1],:) ,'facecolor','interp'); view(2)
axis ([-1 1 -1 1])
% cmin = min(min(v([1:end 1],:)));
% cmax = max(max(v([1:end 1],:)));
% caxis ([cmin cmax]);
title 'sys mean'
 max_v(2) = max(max(v));
min_v(2) = min(min(v));
 caxis ([cmin cmax ]);
 c = colorbar;
%  c.Label.String = 'thickness (mm)';
 %
 
 axis square
 subplot 223
%  shape = shape2 - shape1';
 shape3 = dia_sys_dEPI2ENDOs_plus(:,visualisingMode1)';
colormap jet
shape = shape3(sysPts)';
v = shape(meridian);
surf( X , Y , v([1:end 1],:) ,'facecolor','interp'); view(2)
axis ([-1 1 -1 1])
% cmin = min(min(v([1:end 1],:)));
% cmax = max(max(v([1:end 1],:)));
% caxis ([cmin cmax]);
title ([ 'sys, mode:', num2str(visualisingMode1),' ',num2str(factor1),'std'])
 max_v(3) = max(max(v))
 min_v(3) = min(min(v))

 caxis ([cmin  cmax]);
c = colorbar;
axis square
% c.Label.String = 'thickness (mm)';
%
%
subplot 224
 shape = shape2 - shape1';
%  shape3 = dia_sys_dEPI2ENDOs_plus(:,visualisingMode1)';
%shape = shape3;
colormap jet
shape = shape'
v = shape(meridian);
surf( X , Y , v([1:end 1],:) ,'facecolor','interp'); view(2)
axis ([-1 1 -1 1])
% cmin = min(min(v([1:end 1],:)));
% cmax = max(max(v([1:end 1],:)));
% caxis ([cmin cmax]);
% title ([ 'dia sys, mode:', num2str(visualisingMode1),' ',num2str(factor1),'std'])
title 'difference between negative extreme and mean'
 max_v(3) = max(max(v))
 min_v(3) = min(min(v))

 caxis ([cmin  cmax]);
c = colorbar;
axis square
% c.Label.String = 'thickness (mm)';
%%
