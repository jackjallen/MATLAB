function visualiseModesMovie(data,  dia_sys_myo_mean,principle_dia_sys_myo_eigenvectors, dia_sys_myo_max_b, stage, mode, c, x, y, z)


switch stage
    case 'dia'
        VisualiseDiastolicModes(data, dia_sys_myo_mean, principle_dia_sys_myo_eigenvectors, dia_sys_myo_max_b, stage, mode, c, x, y, z)
    case 'sys'
        VisualiseSystolicModes(data, dia_sys_myo_mean, principle_dia_sys_myo_eigenvectors, dia_sys_myo_max_b, stage, mode, c, x, y, z)
        
end

    function VisualiseDiastolicModes(data, dia_sys_myo_mean, principle_dia_sys_myo_eigenvectors, dia_sys_myo_max_b, stage, mode, c, x, y, z)
        
       figure
        c = c
        
        dia_sys_myo_new_shape(:,mode) = dia_sys_myo_mean(1,:)' + principle_dia_sys_myo_eigenvectors(:,mode)*(c*0.1)*dia_sys_myo_max_b(mode,1);
        
        dia_sys = reshape(dia_sys_myo_new_shape(:,mode), [2*2178 , 3]);
        dia = dia_sys(1:2178,:);
%         
%         x = [-70 ; 40];
%         y = [-65 ; 35];
%         z = [-110 ; 40];
        mode = mode
        
        title (['diastolic, mode ', num2str(mode) ])
        
        xlabel 'x', ylabel 'y', zlabel 'z'
       patch('vertices', dia, 'faces', data(1).diastolic.full.tri, 'facecolor', 'r', 'facealpha', '0.4', 'edgecolor', 'none','FaceLighting','gouraud',  'clipping','off')
        camlight('right')
        %     plot3D(dia)
        view(225, 40)
        axis ([x(1) x(2) y(1) y(2) z(1) z(2)])
        
%         
%         %         subplot 122
%         h2 =  subplot(1,2,2)
%         cla(h2)
%         c = -c
%         
%         dia_sys_myo_new_shape(:,mode) = dia_sys_myo_mean(1,:)' + principle_dia_sys_myo_eigenvectors(:,mode)*c*0.1*dia_sys_myo_max_b(mode,1);
%         
%         dia_sys = reshape(dia_sys_myo_new_shape(:,mode), [2*2178 , 3]);
%         dia = dia_sys(1:2178,:);
% %         
% %         x = [-70 ; 40];
% %         y = [-65 ; 35];
% %         z = [-110 ; 40];
%         mode = mode
%         title (['diastolic, mode ', num2str(mode) ])
%         
%         xlabel 'x', ylabel 'y', zlabel 'z'
%         %      plot3D(dia)
%         patch('vertices', dia, 'faces', data(1).diastolic.full.tri, 'facecolor', 'r', 'facealpha', '0.4', 'edgecolor', 'none','FaceLighting','gouraud')
%         camlight('right')
%         view(225, 40)
%         axis ([x(1) x(2) y(1) y(2) z(1) z(2)])
        
        
    end

    function VisualiseSystolicModes(data, dia_sys_myo_mean, principle_dia_sys_myo_eigenvectors, dia_sys_myo_max_b, stage, mode, c, x, y, z)
        
        c = c
        figure
        dia_sys_myo_new_shape(:,mode) = dia_sys_myo_mean(1,:)' + principle_dia_sys_myo_eigenvectors(:,mode)*c*0.1*dia_sys_myo_max_b(mode,1);
        
        dia_sys = reshape(dia_sys_myo_new_shape(:,mode), [2*2178 , 3]);
        sys = dia_sys(2179:end,:);
        
%         x = [-70 ; 40];
%         y = [-65 ; 35];
%         z = [-110 ; 40];
        mode = mode
        title (['systolic, mode ', num2str(mode) ])
        
        xlabel 'x', ylabel 'y', zlabel 'z'
        patch('vertices', sys, 'faces', data(1).systolic.full.tri, 'facecolor', 'r', 'facealpha', '0.4', 'edgecolor', 'none','FaceLighting','gouraud', 'clipping', 'off') %'specularstrength', 1)
        camlight('right')
     axis tight
        %             plot3D(sys)
        view(225, 40)
    axis ([x(1) x(2) y(1) y(2) z(1) z(2)])
        
%         
%         subplot 122
%         c = -c
%         dia_sys_myo_new_shape(:,mode) = dia_sys_myo_mean(1,:)' + principle_dia_sys_myo_eigenvectors(:,mode)*c*dia_sys_myo_max_b(mode,1);
%         
%         dia_sys = reshape(dia_sys_myo_new_shape(:,mode), [2*2178 , 3]);
%         sys = dia_sys(2179:end,:);
% %         
% %         x = [-70 ; 40];
% %         y = [-65 ; 35];
% %         z = [-110 ; 40];
%         mode = mode
%         title (['systolic, mode ', num2str(mode) ])
%         
%         xlabel 'x', ylabel 'y', zlabel 'z'
%         patch('vertices', sys, 'faces', data(1).systolic.full.tri, 'facecolor', 'r', 'facealpha', '0.4', 'edgecolor', 'none','FaceLighting','gouraud')
%         camlight('right')
%         %       plot3D(sys)
%         view(225, 40)
%         axis ([x(1) x(2) y(1) y(2) z(1) z(2)])
    end







end




