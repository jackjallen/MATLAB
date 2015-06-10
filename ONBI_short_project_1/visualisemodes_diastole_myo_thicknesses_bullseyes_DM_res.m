
diaPts = [1 :1089];
sysPts = [1090: 2178];



R     = linspace( 0 , 1    , size( meridian , 2 )   );
THETA = linspace( 0 , 2*pi , size( meridian , 1 )+1 );
X     = bsxfun( @times , R , cos(THETA).' );
Y     = bsxfun( @times , R , sin(THETA).' );

%%
factor1=nStd
% DETERMINE_dia_sys_dEPI2ENDOs_plus(:,visualisingMode1) = DETERMINE_dia_sys_dEPI2ENDOs_mean(:) + principle_dia_sys_dEPI2ENDOs_eigenvectors(:,visualisingMode1)*max(DETERMINE_dia_sys_dEPI2ENDOs_b(:,visualisingMode1));
% MESA_dia_sys_dEPI2ENDOs_plus(:,visualisingMode1) = MESA_dia_sys_dEPI2ENDOs_mean(:) + principle_dia_sys_dEPI2ENDOs_eigenvectors(:,visualisingMode1)*(factor1)*MESA_dia_sys_dEPI2ENDOs_max_b(visualisingMode1);
factor2=-nStd

%%
figure('name','thickness map')
% for residualIndex = 1:100
residuals(:,visualisingMode1) = MESA_dia_sys_dEPI2ENDOs_max_b(visualisingMode1) - DETERMINE_dia_sys_dEPI2ENDOs_b(:,visualisingMode1); 
residual_map = MESA_dia_sys_dEPI2ENDOs_mean(:) + mean(residuals(:,visualisingMode1))*principle_dia_sys_dEPI2ENDOs_eigenvectors(:,visualisingMode1)
colormap jet
shape = residual_map(diaPts)
% run('make_meridians.m');
v = shape(meridian);
surf( X , Y , v([1:end 1],:) ,'facecolor','interp'); view(2)
axis ([-1 1 -1 1])
title (['from mean( MESA std(mode',num2str(visualisingMode1),') - DETERMINE b values(1:100,mode',num2str(visualisingMode1),') )'])
colorbar
% caxis([ -1 1])
% pause(0.5)
% end
%
% %%
% figure('name', 'PCA on MESA myocardium thickness')
%  subplot 221
% % for CASE = 1:100
% DETERMINE_dia_sys_dEPI2ENDOs_minus(:,visualisingMode1) = DETERMINE_dia_sys_dEPI2ENDOs_mean(:) + principle_dia_sys_dEPI2ENDOs_eigenvectors(:,visualisingMode1)*min(DETERMINE_dia_sys_dEPI2ENDOs_b(:,visualisingMode1));
% MESA_dia_sys_dEPI2ENDOs_minus(:,visualisingMode1) = MESA_dia_sys_dEPI2ENDOs_mean(:) + principle_dia_sys_dEPI2ENDOs_eigenvectors(:,visualisingMode1)*(factor2)*MESA_dia_sys_dEPI2ENDOs_max_b(visualisingMode1);
% 
% visualisingMode1
% 
% clear shape
% 
% shape1 = MESA_dia_sys_dEPI2ENDOs_minus(:,visualisingMode1) -  DETERMINE_dia_sys_dEPI2ENDOs_minus(:,visualisingMode1) ;
% 
% colormap jet
% shape = shape1(diaPts)
% % run('make_meridians.m');
% 
% %%
% v = shape(meridian);
% surf( X , Y , v([1:end 1],:) ,'facecolor','interp'); view(2)
% axis ([-1 1 -1 1])
% % cmin = min(min(v([1:end 1],:)));
% % cmax = max(max(v([1:end 1],:)));
% % caxis ([cmin cmax]);
% title (['diastole, MESA - DETERMINE, mode:', num2str(visualisingMode1),' ',num2str(factor2),'std']);
% max_v(1) = max(max(v));
% min_v(1) = min(min(v));
% %  caxis ([cmin  cmax ]);
%  c = colorbar;
% %  c.Label.String = 'thickness (mm)';
%   axis square
% %   pause
% 
% % end
% subplot 223
% clear shape1
% clear shape
% 
% shape1 = MESA_dia_sys_dEPI2ENDOs_plus(:,visualisingMode1) -  DETERMINE_dia_sys_dEPI2ENDOs_plus(:,visualisingMode1) ;
% 
% colormap jet
% shape = shape1(diaPts)
% % run('make_meridians.m');
% 
% %%
% v = shape(meridian);
% surf( X , Y , v([1:end 1],:) ,'facecolor','interp'); view(2)
% axis ([-1 1 -1 1])
% % cmin = min(min(v([1:end 1],:)));
% % cmax = max(max(v([1:end 1],:)));
% % caxis ([cmin cmax]);
% title (['diastole, MESA - DETERMINE, mode:', num2str(visualisingMode1),' ',num2str(factor1),'std']);
% max_v(1) = max(max(v));
% min_v(1) = min(min(v));
% %  caxis ([cmin  cmax ]);
%  c = colorbar;
% %  c.Label.String = 'thickness (mm)';
%   axis square
%   
%   
% %
% %  shape2 =dia_sys_dEPI2ENDOs_mean';
% % % figure('name', 'PCA on myocardium thicknesses')
% %  subplot 222
% % colormap jet
% % shape = shape2(diaPts)
% % % run('make_meridians.m');
% % v = shape(meridian);
% % surf( X , Y , v([1:end 1],:) ,'facecolor','interp'); view(2)
% % axis ([-1 1 -1 1])
% % % cmin = min(min(v([1:end 1],:)));
% % % cmax = max(max(v([1:end 1],:)));
% % % caxis ([cmin cmax]);
% % title 'dia mean'
% %  max_v(2) = max(max(v));
% % min_v(2) = min(min(v));
% %  caxis ([cmin cmax ]);
% %  c = colorbar;
% % %  c.Label.String = 'thickness (mm)';
% %  %
% %  
% %  axis square
% %  subplot 223
% % %  shape = shape2 - shape1';
% %  shape3 = dia_sys_dEPI2ENDOs_plus(:,visualisingMode1)';
% % colormap jet
% % shape = shape3(diaPts);
% % v = shape(meridian);
% % surf( X , Y , v([1:end 1],:) ,'facecolor','interp'); view(2)
% % axis ([-1 1 -1 1])
% % % cmin = min(min(v([1:end 1],:)));
% % % cmax = max(max(v([1:end 1],:)));
% % % caxis ([cmin cmax]);
% % title ([ 'dia , mode:', num2str(visualisingMode1),' ',num2str(factor1),'std'])
% %  max_v(3) = max(max(v))
% %  min_v(3) = min(min(v))
% % 
% %  caxis ([cmin  cmax]);
% % c = colorbar;
% % axis square
% % % c.Label.String = 'thickness (mm)';
% % %
% % %
% % subplot 224
% %  shape = shape2 - shape1';
% % %  shape3 = dia_sys_dEPI2ENDOs_plus([diaPts;sysPts],visualisingMode1)';
% % %shape = shape3;
% % colormap jet
% % shape = shape(diaPts)'
% % v = shape(meridian);
% % surf( X , Y , v([1:end 1],:) ,'facecolor','interp'); view(2)
% % axis ([-1 1 -1 1])
% % % cmin = min(min(v([1:end 1],:)));
% % % cmax = max(max(v([1:end 1],:)));
% % % caxis ([cmin cmax]);
% % % title ([ 'dia sys, mode:', num2str(visualisingMode1),' ',num2str(factor1),'std'])
% % title 'difference between negative extreme and mean'
% % 
% %  max_v(3) = max(max(v))
% %  min_v(3) = min(min(v))
% % 
% %  caxis ([cmin  cmax]);
% % c = colorbar;
% % axis square
% % % c.Label.String = 'thickness (mm)';
% % 
