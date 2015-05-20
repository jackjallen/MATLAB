function plotVolumes(data, sys_epi_volumes, sys_endo_volumes, dia_epi_volumes, dia_endo_volumes)

%% plot volume histograms - comparing DETERMINE with MESA
% 1mm^3 = 0.001ml
%diastolic endocardium volumes


for d =  data(1).DETERMINE_indices'
    
    DETERMINE_dia_endo_volumes(d) = dia_endo_volumes(d);
    DETERMINE_dia_epi_volumes(d) = dia_epi_volumes(d);
    DETERMINE_sys_endo_volumes(d) = sys_endo_volumes(d);
    DETERMINE_sys_epi_volumes(d) = sys_epi_volumes(d);
end
for m = data(1).MESA_indices'
    
    MESA_dia_endo_volumes(m) = dia_endo_volumes(m);
    MESA_dia_epi_volumes(m) = dia_epi_volumes(m);
    MESA_sys_endo_volumes(m) = sys_endo_volumes(m);
    MESA_sys_epi_volumes(m) = sys_epi_volumes(m);
end


DETERMINE_dia_endo_volumes( :,all(~DETERMINE_dia_endo_volumes,1) ) = [];
DETERMINE_dia_epi_volumes( :,all(~DETERMINE_dia_epi_volumes,1) ) = [];
DETERMINE_sys_endo_volumes( :,all(~DETERMINE_sys_endo_volumes,1) ) = [];
DETERMINE_sys_epi_volumes( :,all(~DETERMINE_sys_epi_volumes,1) ) = [];

MESA_dia_endo_volumes( :,all(~MESA_dia_endo_volumes,1) ) = [];
MESA_dia_epi_volumes( :,all(~MESA_dia_epi_volumes,1) ) = [];
MESA_sys_endo_volumes( :,all(~MESA_sys_endo_volumes,1) ) = [];
MESA_sys_epi_volumes( :,all(~MESA_sys_epi_volumes,1) ) = [];



figure
subplot 121
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

subplot 122
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

end