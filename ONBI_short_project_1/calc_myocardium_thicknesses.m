% compute myocardium wall thickness and basic statistics

for i = 1:401
    %  wall thicknesses
    [~,~,data(i).dia_dEPI2ENDO] = vtkClosestElement( data(i).diastolic.endo , data(i).diastolic.epi.xyz );
    [~,~,data(i).sys_dEPI2ENDO] = vtkClosestElement( data(i).systolic.endo , data(i).systolic.epi.xyz );
    
    diastolicEndos(i,:) = data(i).diastolic.endo.xyz(:);
    systolicEndos(i,:)= data(i).systolic.endo.xyz(:);
    diastolicEpis(i,:) = data(i).diastolic.epi.xyz(:);
    systolicEpis(i,:)= data(i).systolic.epi.xyz(:);
    
    dia_dEPI2ENDOs(i,:) = data(i).dia_dEPI2ENDO;
    sys_dEPI2ENDOs(i,:) = data(i).sys_dEPI2ENDO;
 
    data(i).dia_dEPI2ENDO_vars = var(data(i).dia_dEPI2ENDO);
    data(i).sys_dEPI2ENDO_vars = var(data(i).sys_dEPI2ENDO);
    
    data(i).dia_sys_dEPI2ENDO_vars = var([data(i).dia_dEPI2ENDO ; data(i).sys_dEPI2ENDO ])
    
    data(i).dia_dEPI2ENDO_stds = sqrt(data(i).dia_dEPI2ENDO_vars);
    data(i).sys_dEPI2ENDO_stds = sqrt(data(i).sys_dEPI2ENDO_vars);
     
    data(i).myo_T_changes = data(i).sys_dEPI2ENDO - data(i).dia_dEPI2ENDO;
   
    diaMeanT(i) = mean(data(i).dia_dEPI2ENDO);
    sysMeanT(i) = mean(data(i).sys_dEPI2ENDO);
    
    diaMedianT(i) = median(data(i).dia_dEPI2ENDO);
    sysMedianT(i) = median(data(i).sys_dEPI2ENDO);
    
    diaModeT(i) = mode(data(i).dia_dEPI2ENDO);
    sysModeT(i) = mode(data(i).sys_dEPI2ENDO);
    
    deltaMeanT(i) = mean(data(i).sys_dEPI2ENDO) - mean(data(i).dia_dEPI2ENDO);
    
end

DETERMINEmeanDiaEndo = mean(diastolicEndos(data(1).DETERMINE_indices,:),1);
DETERMINEmeanDiaEpi = mean(diastolicEpis(data(1).DETERMINE_indices,:),1);
DETERMINEmeanDiaT = mean(dia_dEPI2ENDOs(data(1).DETERMINE_indices,:),1);
MESAmeanDiaEndo = mean(diastolicEndos(data(1).MESA_indices,:),1);
MESAmeanDiaEpi = mean(diastolicEpis(data(1).MESA_indices,:),1);
MESAmeanDiaT = mean(dia_dEPI2ENDOs(data(1).MESA_indices,:),1);

DETERMINEmeanSysEndo = mean(systolicEndos(data(1).DETERMINE_indices,:),1);
DETERMINEmeanSysEpi = mean(systolicEpis(data(1).DETERMINE_indices,:),1);
DETERMINEmeanSysT = mean(sys_dEPI2ENDOs(data(1).DETERMINE_indices,:),1);
MESAmeanSysEndo = mean(systolicEndos(data(1).MESA_indices,:),1);
MESAmeanSysEpi = mean(systolicEpis(data(1).MESA_indices,:),1);
MESAmeanSysT = mean(sys_dEPI2ENDOs(data(1).MESA_indices,:),1);

DETERMINEdiaMyoThicknessVariance = cell2mat({data(data(1).DETERMINE_indices(:)).dia_dEPI2ENDO_vars}) ;
MESAdiaMyoThicknessVariance = cell2mat({data(data(1).MESA_indices).dia_dEPI2ENDO_vars});
DETERMINEsysMyoThicknessVariance = cell2mat({data(data(1).DETERMINE_indices).sys_dEPI2ENDO_vars});
MESAsysMyoThicknessVariance = cell2mat({data(data(1).MESA_indices).sys_dEPI2ENDO_vars});

DETERMINEmyoTchangesVar = var(cell2mat({data(data(1).DETERMINE_indices).myo_T_changes}))
MESAmyoTchangesVar = var(cell2mat({data(data(1).MESA_indices).myo_T_changes}))

DETERMINEdeltaMeanT = deltaMeanT(data(1).DETERMINE_indices)
MESAdeltaMeanT = deltaMeanT(data(1).MESA_indices)

DETERMINE_DiaSys_dEPI2ENDO_vars = cell2mat({data(data(1).DETERMINE_indices).dia_sys_dEPI2ENDO_vars})
MESA_DiaSys_dEPI2ENDO_vars = cell2mat({data(data(1).MESA_indices).dia_sys_dEPI2ENDO_vars})