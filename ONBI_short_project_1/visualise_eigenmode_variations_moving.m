%% visualise mode variations
visualisingMode = 2
stage = 'sys'
%animation
for n= 1
    
    x = [-60 ; 50];
    y = [-55 ; 50];
    z = [-110 ; 15];
    loops = 10;
    F(loops) = struct('cdata',[],'colormap',[]);
    
    for f = 1:10
        tmpf = sort(1:10,'descend')
        c = tmpf(f); % -1, -0.9, ... , -0.1
        pause(0.1)
        
        visualiseModesMovie(data, dia_sys_myo_mean,principle_dia_sys_myo_eigenvectors, dia_sys_myo_max_b, stage , visualisingMode, c,x,y,z)
        drawnow
        F(f) = getframe;
    end
    for f = 1:10
        
        c = -f % 0.1, 0.2, ... , 1
        pause(0.1)
        hold off
        visualiseModesMovie(data, dia_sys_myo_mean,principle_dia_sys_myo_eigenvectors, dia_sys_myo_max_b, stage, visualisingMode, c,x,y,z)
        drawnow
        F(f+10) = getframe;
    end
end
