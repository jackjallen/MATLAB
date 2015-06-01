subplot 332
factor=c
dia_sys_myo_new_shape_plus(:,visualisingMode2) = dia_sys_myo_mean(:) + principle_dia_sys_myo_eigenvectors(:,visualisingMode2)*(factor)*dia_sys_myo_max_b(visualisingMode2,1);
dia_sys = reshape(dia_sys_myo_new_shape_plus(:,visualisingMode2), [2*2178 , 3]);
dia = dia_sys(1:2178,:);
%         x = [-70 ; 40];
%         y = [-65 ; 35];
%         z = [-110 ; 40];
title (['Mode ', num2str(visualisingMode2),', ',num2str(factor),'std' ])
% tri = data(1).diastolic.full.tri
xlabel 'x', ylabel 'y', zlabel 'z'
patch('vertices', dia, 'faces', data(1).diastolic.full.tri, 'facecolor', 'r', 'facealpha', '0.4', 'edgecolor', 'none','FaceLighting','gouraud', 'clipping','off')
camlight('right')
%     plot3D(dia)
view(225, 40)
axis ([x(1) x(2) y(1) y(2) z(1) z(2)])

subplot 334
factor = -c
dia_sys_myo_new_shape_minus(:,visualisingMode1) = dia_sys_myo_mean(:) + principle_dia_sys_myo_eigenvectors(:,visualisingMode1)*(factor)*dia_sys_myo_max_b(visualisingMode1,1);
shape = reshape(dia_sys_myo_new_shape_minus(:,visualisingMode1), [2*2178 , 3]);
shape = shape(1:2178,:)
patch('vertices', shape, 'faces', data(1).diastolic.full.tri, 'facecolor', 'r', 'facealpha', '0.4', 'edgecolor', 'none','FaceLighting','gouraud', 'clipping','off')
camlight('right')
%     plot3D(dia)
title (['Mode ', num2str(visualisingMode1),', ',num2str(factor),'std' ])
view(225, 40)
axis ([x(1) x(2) y(1) y(2) z(1) z(2)])

subplot 335
shape = reshape(dia_sys_myo_mean(:), [2*2178 , 3]);
shape = shape(1:2178,:)
patch('vertices', shape, 'faces', data(1).diastolic.full.tri, 'facecolor', 'r', 'facealpha', '0.4', 'edgecolor', 'none','FaceLighting','gouraud', 'clipping','off')
camlight('right')
%     plot3D(dia)
title (['mean' ])
view(225, 40)
axis ([x(1) x(2) y(1) y(2) z(1) z(2)])



subplot 336
factor = c
dia_sys_myo_new_shape_plus(:,visualisingMode1) = dia_sys_myo_mean(:) + principle_dia_sys_myo_eigenvectors(:,visualisingMode1)*factor*dia_sys_myo_max_b(visualisingMode1,1);
shape = reshape(dia_sys_myo_new_shape_plus(:,visualisingMode1), [2*2178 , 3]);
dia = shape(1:2178,:);
%         x = [-70 ; 40];
%         y = [-65 ; 35];
%         z = [-110 ; 40];
title (['Mode ', num2str(visualisingMode1),', ',num2str(factor),'std' ])
xlabel 'x', ylabel 'y', zlabel 'z'
%      plot3D(dia)
patch('vertices', dia, 'faces', data(1).diastolic.full.tri, 'facecolor', 'r', 'facealpha', '0.4', 'edgecolor', 'none','FaceLighting','gouraud',  'clipping','off')
camlight('right')
view(225, 40)
axis ([x(1) x(2) y(1) y(2) z(1) z(2)])

subplot 338
factor = -c
dia_sys_myo_new_shape_minus(:,visualisingMode2) = dia_sys_myo_mean(:) + principle_dia_sys_myo_eigenvectors(:,visualisingMode2)*factor*dia_sys_myo_max_b(visualisingMode2,1);
shape = reshape(dia_sys_myo_new_shape_minus(:,visualisingMode2), [2*2178 , 3]);
dia = shape(1:2178,:);
%         x = [-70 ; 40];
%         y = [-65 ; 35];
%         z = [-110 ; 40];
title (['Mode ', num2str(visualisingMode2),', ',num2str(factor),'std' ])
xlabel 'x', ylabel 'y', zlabel 'z'
%      plot3D(dia)
patch('vertices', dia, 'faces', data(1).diastolic.full.tri, 'facecolor', 'r', 'facealpha', '0.4', 'edgecolor', 'none','FaceLighting','gouraud',  'clipping','off')
camlight('right')
view(225, 40)
axis ([x(1) x(2) y(1) y(2) z(1) z(2)])