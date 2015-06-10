%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot 331
factor1 = -nStd;
factor2 = nStd;
dia_sys_dEPI2ENDOs_new_shape_minus_plus(:,visualisingMode1) = dia_sys_dEPI2ENDOs_mean(:) + principle_dia_sys_dEPI2ENDOs_eigenvectors(:,visualisingMode1)*(factor1)*dia_sys_dEPI2ENDOs_max_b(visualisingMode1,1) + ...
principle_dia_sys_dEPI2ENDOs_eigenvectors(:,visualisingMode2)*(factor2)*dia_sys_dEPI2ENDOs_max_b(visualisingMode2,1) ;
shape = reshape(dia_sys_dEPI2ENDOs_new_shape_minus_plus(:,visualisingMode1),[4356 3]);
shape = shape(2179:end,:);
patch('vertices', shape, 'faces', data(1).diastolic.full.tri, 'facecolor', 'r', 'facealpha', '0.4', 'edgecolor', 'none','FaceLighting','gouraud', 'clipping','off')
camlight('right')
%     plot3D(dia)
title (['Mode ', num2str(visualisingMode1),', ',num2str(factor1),'std',', ','Mode ', num2str(visualisingMode2),', ',num2str(factor2),'std' ])
view(225, 40)
axis ([x(1) x(2) y(1) y(2) z(1) z(2)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot 333
factor1 = nStd;
factor2 = nStd;
dia_sys_dEPI2ENDOs_new_shape_plus_plus(:,visualisingMode1) = dia_sys_dEPI2ENDOs_mean(:) + principle_dia_sys_dEPI2ENDOs_eigenvectors(:,visualisingMode1)*(factor1)*dia_sys_dEPI2ENDOs_max_b(visualisingMode1,1) + ...
principle_dia_sys_dEPI2ENDOs_eigenvectors(:,visualisingMode2)*(factor2)*dia_sys_dEPI2ENDOs_max_b(visualisingMode2,1) ;
shape = reshape(dia_sys_dEPI2ENDOs_new_shape_plus_plus(:,visualisingMode1), [2*2178 , 3]);
shape = shape(2179:end,:);
%
%         x = [-70 ; 40];
%         y = [-65 ; 35];
%         z = [-110 ; 40];
title (['Mode ', num2str(visualisingMode1),', ',num2str(factor1),'std',', ','Mode ', num2str(visualisingMode2),', ',num2str(factor2),'std' ])
% tri = data(1).diastolic.full.tri
xlabel 'x', ylabel 'y', zlabel 'z'
patch('vertices', shape, 'faces', data(1).diastolic.full.tri, 'facecolor', 'r', 'facealpha', '0.4', 'edgecolor', 'none','FaceLighting','gouraud', 'clipping','off')
camlight('right')
%     plot3D(dia)
view(225, 40)
axis ([x(1) x(2) y(1) y(2) z(1) z(2)])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot 337
factor1 = -nStd;
factor2 = -nStd;
dia_sys_dEPI2ENDOs_new_shape_minus_minus(:,visualisingMode1) = dia_sys_dEPI2ENDOs_mean(:) + principle_dia_sys_dEPI2ENDOs_eigenvectors(:,visualisingMode1)*(factor1)*dia_sys_dEPI2ENDOs_max_b(visualisingMode1,1) + ...
principle_dia_sys_dEPI2ENDOs_eigenvectors(:,visualisingMode2)*(factor2)*dia_sys_dEPI2ENDOs_max_b(visualisingMode2,1) ;
shape = reshape(dia_sys_dEPI2ENDOs_new_shape_minus_minus(:,visualisingMode1), [2*2178 , 3]);
shape = shape(2179:end,:);
%         x = [-70 ; 40];
%         y = [-65 ; 35];
%         z = [-110 ; 40];
title (['Mode ', num2str(visualisingMode1),', ',num2str(factor1),'std',', ','Mode ', num2str(visualisingMode2),', ',num2str(factor2),'std' ])
xlabel 'x', ylabel 'y', zlabel 'z'
%      plot3D(dia)
patch('vertices', shape, 'faces', data(1).diastolic.full.tri, 'facecolor', 'r', 'facealpha', '0.4', 'edgecolor', 'none','FaceLighting','gouraud',  'clipping','off')
camlight('right')
view(225, 40)
axis ([x(1) x(2) y(1) y(2) z(1) z(2)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot 339
factor1 = nStd;
factor2 = -nStd;
dia_sys_dEPI2ENDOs_new_shape_plus_minus(:,visualisingMode1) = dia_sys_dEPI2ENDOs_mean(:) + principle_dia_sys_dEPI2ENDOs_eigenvectors(:,visualisingMode1)*(factor1)*dia_sys_dEPI2ENDOs_max_b(visualisingMode1,1) + ...
principle_dia_sys_dEPI2ENDOs_eigenvectors(:,visualisingMode2)*(factor2)*dia_sys_dEPI2ENDOs_max_b(visualisingMode2,1) ;
shape = reshape(dia_sys_dEPI2ENDOs_new_shape_plus_minus(:,visualisingMode1), [2*2178 , 3]);
shape = shape(2179:end,:);
%         x = [-70 ; 40];
%         y = [-65 ; 35];
%         z = [-110 ; 40];
title (['Mode ', num2str(visualisingMode1),', ',num2str(factor1),'std',', ','Mode ', num2str(visualisingMode2),', ',num2str(factor2),'std' ])
xlabel 'x', ylabel 'y', zlabel 'z'
%      plot3D(dia)
patch('vertices', shape, 'faces', data(1).diastolic.full.tri, 'facecolor', 'r', 'facealpha', '0.4', 'edgecolor', 'none','FaceLighting','gouraud',  'clipping','off')
camlight('right')
view(225, 40)
axis ([x(1) x(2) y(1) y(2) z(1) z(2)])