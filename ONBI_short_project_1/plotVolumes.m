function plotVolumes(data, sys_epi_volumes, sys_endo_volumes, dia_epi_volumes, dia_endo_volumes)

%% plot volume histograms - comparing DETERMINE with MESA
% 1mm^3 = 0.001ml
%diastolic endocardium volumes

for i = 1:400
    for d = 1:100
        if data(1).DETERMINE_indices(d)==i
            
            DETERMINE_dia_endo_volumes(d) = dia_endo_volumes(i);
            DETERMINE_dia_epi_volumes(d) = dia_epi_volumes(i);
            DETERMINE_sys_endo_volumes(d) = sys_endo_volumes(i);
            DETERMINE_sys_epi_volumes(d) = sys_epi_volumes(i);
        end
        if data(1).MESA_indices(d)==i
            
            MESA_dia_endo_volumes(d) = dia_endo_volumes(i);
            MESA_dia_epi_volumes(d) = dia_epi_volumes(i);
            MESA_sys_endo_volumes(d) = sys_endo_volumes(i);
            MESA_sys_epi_volumes(d) = sys_epi_volumes(i);
        end
        
    end
    
end
figure
hold on
nbins = 25;
histogram(DETERMINE_dia_endo_volumes*0.001,nbins)
histogram(MESA_dia_endo_volumes*0.001,nbins)
legend 'DETERMINE' ' MESA'
title 'diastolic endocardium volumes'
xlabel 'endocardium volume (ml)'
ylabel 'frequency'
% print('compare diastolic endocardium volumes','-dpng')

% systolic endocardium volumes
figure
hold on
nbins = 25;
histogram(DETERMINE_sys_endo_volumes*0.001,nbins)
histogram(MESA_sys_endo_volumes*0.001,nbins)
legend 'DETERMINE' ' MESA'
title 'systolic endocardium volumes'
xlabel 'endocardium volume (ml)'
ylabel 'frequency'
% print('compare systolic endocardium volumes','-dpng')

%calculate endocardium ejection fractions EF
% SV = diastolic volume - systolic volume
% EF = (SV/(diastolic volume))*100