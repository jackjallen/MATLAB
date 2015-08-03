trainingdata(:,1) = [DETERMINE_ejectionFractions' ; MESA_ejectionFractions'];
%modes of variation
trainingdata(:,2:6) = [DETERMINE_sys_myo_b(:,1:5) ; MESA_sys_myo_b(:,1:5)]; %from the PDM
trainingdata(:,7:11) = [DETERMINE_dia_myo_b(:,1:5) ; MESA_dia_myo_b(:,1:5)];
trainingdata(:,12:26) = [DETERMINE_dia_sys_myo_b(:,1:15) ; MESA_dia_sys_myo_b(:,1:15)];
%sphericity measure
trainingdata(:,27) = [DETERMINE_dia_endo_AVratios ; MESA_dia_endo_AVratios ];
trainingdata(:,28) = [DETERMINE_sys_endo_AVratios ; MESA_sys_endo_AVratios ];
trainingdata(:,29) = [DETERMINE_dia_epi_AVratios ; MESA_dia_epi_AVratios ];
trainingdata(:,30) = [DETERMINE_sys_epi_AVratios ; MESA_sys_epi_AVratios ];
% mean thickness
trainingdata(:,31) = [diaMeanT(data(1).DETERMINE_indices)' ; diaMeanT(data(1).MESA_indices)'];
trainingdata(:,32) = [sysMeanT(data(1).DETERMINE_indices)' ; sysMeanT(data(1).MESA_indices)'];
% mean(sys thickness) - mean(dia thickness) e.g. change in mean thickness from diastole to systole
trainingdata(:,33) = [deltaMeanT(data(1).DETERMINE_indices)' ; deltaMeanT(data(1).MESA_indices)'];
%most common thickness
trainingdata(:,34) = [diaModeT(data(1).DETERMINE_indices)' ; diaModeT(data(1).MESA_indices)'];
trainingdata(:,35) = [sysModeT(data(1).DETERMINE_indices)' ; sysModeT(data(1).MESA_indices)'];
% median thickness
trainingdata(:,36) = [diaMedianT(data(1).DETERMINE_indices)' ; diaMedianT(data(1).MESA_indices)'];
trainingdata(:,37) = [sysMedianT(data(1).DETERMINE_indices)' ; sysMedianT(data(1).MESA_indices)'];
%thickness variance for diastole and systole seperately
trainingdata(:,38) = [DETERMINEdiaMyoThicknessVariance' ; MESAdiaMyoThicknessVariance' ];
trainingdata(:,39) = [DETERMINEsysMyoThicknessVariance' ; MESAsysMyoThicknessVariance' ];
% var(sys thickness - dia thickness)  e.g. variance in the difference between diastole and systole
trainingdata(:,40) = [DETERMINEmyoTchangesVar' ;MESAmyoTchangesVar'];
% var([ dia thickness ; sys thickness])
 trainingdata(:,41) = [DETERMINE_DiaSys_dEPI2ENDO_vars' ; MESA_DiaSys_dEPI2ENDO_vars' ];
%variance of the thickness in each of the 3 different sections
trainingdata(:,42) = [cell2mat({data(data(1).DETERMINE_indices).dia_vA_var})'  ; cell2mat({data(data(1).MESA_indices).dia_vA_var})'];
trainingdata(:,43) = [cell2mat({data(data(1).DETERMINE_indices).dia_vB_var})'  ; cell2mat({data(data(1).MESA_indices).dia_vB_var})'];
trainingdata(:,44) = [cell2mat({data(data(1).DETERMINE_indices).dia_vC_var})'  ; cell2mat({data(data(1).MESA_indices).dia_vC_var})'];
trainingdata(:,45) = [cell2mat({data(data(1).DETERMINE_indices).sys_vA_var})'  ; cell2mat({data(data(1).MESA_indices).sys_vA_var})'];
trainingdata(:,46) = [cell2mat({data(data(1).DETERMINE_indices).sys_vB_var})'  ; cell2mat({data(data(1).MESA_indices).sys_vB_var})'];
trainingdata(:,47) = [cell2mat({data(data(1).DETERMINE_indices).sys_vC_var})'  ; cell2mat({data(data(1).MESA_indices).sys_vC_var})'];
trainingdata(:,48) = [cell2mat({data(data(1).DETERMINE_indices).dia_sys_vA_var})'  ; cell2mat({data(data(1).MESA_indices).dia_sys_vA_var})'];
trainingdata(:,49) = [cell2mat({data(data(1).DETERMINE_indices).dia_sys_vB_var})'  ; cell2mat({data(data(1).MESA_indices).dia_sys_vB_var})'];
trainingdata(:,50) = [cell2mat({data(data(1).DETERMINE_indices).dia_sys_vC_var})'  ; cell2mat({data(data(1).MESA_indices).dia_sys_vC_var})'];
%mean thickness for each of the 3 different sections
trainingdata(:,51) = [cell2mat({data(data(1).DETERMINE_indices).dia_vA_mean})'  ; cell2mat({data(data(1).MESA_indices).dia_vA_mean})'];
trainingdata(:,52) = [cell2mat({data(data(1).DETERMINE_indices).dia_vB_mean})'  ; cell2mat({data(data(1).MESA_indices).dia_vB_mean})'];
trainingdata(:,53) = [cell2mat({data(data(1).DETERMINE_indices).dia_vC_mean})'  ; cell2mat({data(data(1).MESA_indices).dia_vC_mean})'];
trainingdata(:,54) = [cell2mat({data(data(1).DETERMINE_indices).sys_vA_mean})'  ; cell2mat({data(data(1).MESA_indices).sys_vA_mean})'];
trainingdata(:,55) = [cell2mat({data(data(1).DETERMINE_indices).sys_vB_mean})'  ; cell2mat({data(data(1).MESA_indices).sys_vB_mean})'];
trainingdata(:,56) = [cell2mat({data(data(1).DETERMINE_indices).sys_vC_mean})'  ; cell2mat({data(data(1).MESA_indices).sys_vC_mean})'];
trainingdata(:,57) = [cell2mat({data(data(1).DETERMINE_indices).dia_sys_vA_mean})'  ; cell2mat({data(data(1).MESA_indices).dia_sys_vA_mean})'];
trainingdata(:,58) = [cell2mat({data(data(1).DETERMINE_indices).dia_sys_vB_mean})'  ; cell2mat({data(data(1).MESA_indices).dia_sys_vB_mean})'];
trainingdata(:,59) = [cell2mat({data(data(1).DETERMINE_indices).dia_sys_vC_mean})'  ; cell2mat({data(data(1).MESA_indices).dia_sys_vC_mean})'];
trainingdata(:,60) = [sys_endo_volumes(data(1).DETERMINE_indices) ; sys_endo_volumes(data(1).MESA_indices) ];
trainingdata(:,61) = [sys_epi_volumes(data(1).DETERMINE_indices) ; sys_epi_volumes(data(1).MESA_indices) ];
trainingdata(:,62) = [dia_endo_volumes(data(1).DETERMINE_indices) ; dia_endo_volumes(data(1).MESA_indices) ];
trainingdata(:,63) = [dia_epi_volumes(data(1).DETERMINE_indices) ; dia_epi_volumes(data(1).MESA_indices) ];
trainingdata(:,64) = [log(sys_endo_volumes(data(1).DETERMINE_indices)) ; log(sys_endo_volumes(data(1).MESA_indices)) ];
trainingdata(:,65) = [log(sys_epi_volumes(data(1).DETERMINE_indices)) ; log(sys_epi_volumes(data(1).MESA_indices)) ];
trainingdata(:,66) = [log(dia_endo_volumes(data(1).DETERMINE_indices)) ; log(dia_endo_volumes(data(1).MESA_indices)) ];
trainingdata(:,67) = [log(dia_epi_volumes(data(1).DETERMINE_indices)) ; log(dia_epi_volumes(data(1).MESA_indices)) ];
% sum(systole thickness)/sum(diastole thickness
trainingdata(:,68) = [MyoThicknessSysDiaRatio(data(1).DETERMINE_indices)' ; MyoThicknessSysDiaRatio(data(1).MESA_indices)'];


%% training data labels
trainingdataLabels = [
    'EF                ' ;...
    'sys myo b1        ' ;...
    'sys myo b2        ' ;...
    'sys myo b3        ' ;...
    'sys myo b4        ' ;...
    'sys myo b5        ' ;...
    'dia myo b1        ' ;...
    'dia myo b2        ' ;...
    'dia myo b3        ' ;...
    'dia myo b4        ' ;...
    'dia myo b5        ' ;...
    'dia sys myo b1    ' ;...
    'dia sys myo b2    ' ;...
    'dia sys myo b3    ' ;...
    'dia sys myo b4    ' ;...
    'dia sys myo b5    ' ;...
    'dia sys myo b6    ' ;...
    'dia sys myo b7    ' ;...
    'dia sys myo b8    ' ;...
    'dia sys myo b9    ' ;...
    'dia sys myo b10   ' ;...
    'dia sys myo b11   ' ;...
    'dia sys myo b12   ' ;...
    'dia sys myo b13   ' ;...
    'dia sys myo b14   ' ;...
    'dia sys myo b15   ' ;...
    'diaendo sphericity' ;...
    'sysendo sphericity' ;...
    'diaepi sphericity ' ;...
    'sysepi sphericity ' ;...
    'diaMeanT          ' ;...
    'sysMeanT          ' ;...
    'deltameanT        ' ;...
    'diaModeT          ' ;...
    'sysModeT          ' ;...
    'diaMedianT        ' ;...
    'sysMedianT        ' ;...
    'diaMyoThicknessVar' ;...
    'sysMyoThicknessVar' ;...
    'myoTchangesVar    ' ;...
    'DiaSysdEPI2ENDOvar' ;...
    'dia_vA_var        ' ;...
    'dia_vB_var        ' ;...
    'dia_vC_var        ' ;...
    'sys_vA_var        ' ;...
    'sys_vB_var        ' ;...
    'sys_vC_var        ' ;...
    'dia_sys_vA_var    ' ;...
    'dia_sys_vB_var    ' ;...
    'dia_sys_vC_var    ' ;...
    'dia_vA_mean       ' ;...
    'dia_vB_mean       ' ;...
    'dia_vC_mean       ' ;...
    'sys_vA_mean       ' ;...
    'sys_vB_mean       ' ;...
    'sys_vC_mean       ' ;...
    'dia_sys_vA_mean   ' ;...
    'dia_sys_vB_mean   ' ;...
    'dia_sys_vC_mean   ' ;...
    'sys endo volume   ' ;...
    'sys epi volume    ' ;...
    'dia endo volume   ' ;...%62
    'dia epi volume    ' ;...%63
    'log sys endo vol  ' ;...%64
    'log sys epi vol   ' ;...%65
    'log dia endo vol  ' ;...%66
    'log dia epi vol   ' ;...%67
    'sum(sysT)/sum(diaT' ];  %68
    
    
    

%% class names
names = char(200,1);
names(1:100,1) = 'd';
names(101:200,1) = 'm';