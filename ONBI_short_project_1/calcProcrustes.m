function [data, dia_myo_mean, sys_myo_mean, dia_sys_myo_mean, all_training_diastolic_myo_shapes, all_training_systolic_myo_shapes , all_training_dia_sys_myo_shapes] = calcProcrustes(data, pI, dia_myo_reference, sys_myo_reference, dia_sys_myo_reference)

all_training_diastolic_myo_shapes = zeros(401, 6534);
all_training_systolic_myo_shapes = zeros(401, 6534);
all_training_dia_sys_myo_shapes = zeros(401, 13068);

concatIndices = sort([data(1).MESA_indices ; data(1).DETERMINE_indices ]);

for p = 1:pI
    
    % iterate so that procrustes is performed on each case (400 patients).
    
    %     for n = 1:200
    for   i = concatIndices(1:end-1)';
        % % transform endo and epi (individually)
        % [data(i).diastolic.endo.procrustes_d, data(i).diastolic.endo.xyz] = procrustes(dia_endo_reference, data(i).diastolic.endo.xyz(:));
        % [data(i).diastolic.epi.procrustes_d, data(i).diastolic.epi.xyz] = procrustes(dia_epi_reference, data(i).diastolic.epi.xyz(:));
        % [data(i).systolic.endo.procrustes_d, data(i).systolic.endo.xyz] = procrustes(sys_endo_reference, data(i).systolic.endo.xyz(:));
        % [data(i).systolic.epi.procrustes_d, data(i).systolic.epi.xyz] = procrustes(sys_epi_reference, data(i).systolic.epi.xyz(:));
        
        % transform endo and epi as one shape vector
        [data(i).diastolic.myo.procrustes_d, diastolic_myo_shapes(i).xyz, dia_myo_transform(i)] = procrustes(dia_myo_reference, data(i).diastolic.myo.xyz,'scaling',false,'reflection',false);
        [data(i).systolic.myo.procrustes_d, systolic_myo_shapes(i).xyz, sys_myo_transform(i)] = procrustes(sys_myo_reference, data(i).systolic.myo.xyz,'scaling',false,'reflection',false); %,'scaling',false); % with scaling SMM0362 goes much smaller than the others...
        
        [data(i).dia_sys.myo.procrustes_d, dia_sys_myo_shapes(i).xyz, dia_sys_transform(i)] = procrustes(dia_sys_myo_reference, data(i).dia_sys.myo.xyz, 'scaling', false, 'reflection', false);
        
        all_training_diastolic_myo_shapes(i,:) = diastolic_myo_shapes(i).xyz(:);
        all_training_systolic_myo_shapes(i,:) = systolic_myo_shapes(i).xyz(:);
        
        all_training_dia_sys_myo_shapes(i,:) = dia_sys_myo_shapes(i).xyz(:);
        
        %         end
        %sums for finding new means later on
        % dia_endo_sum = data(i).diastolic.endo.xyz + dia_endo_sum;
        % dia_epi_sum = data(i).diastolic.epi.xyz + dia_epi_sum;
        % sys_endo_sum = data(i).systolic.endo.xyz + sys_endo_sum;
        % sys_epi_sum = data(i).systolic.epi.xyz + sys_epi_sum;
        
        % dia_myo_sum = data(i).diastolic.myo.xyz + dia_myo_sum;
        % sys_myo_sum = data(i).systolic.myo.xyz + sys_myo_sum;
        %      end
        
        
        
        %         for i = data(1).MESA_indices'
        %
        %             MESA_diastolic_myo_shapes(n,:) = all_training_diastolic_myo_shapes(i,:);
        %             MESA_systolic_myo_shapes(n,:) = all_training_systolic_myo_shapes(i,:);
        %
        %             MESA_dia_sys_myo_shapes(n,:) =  all_training_dia_sys_myo_shapes(i,:);
        %
        %         end
        %         for i = data(1).DETERMINE_indices'
        %
        %             DETERMINE_diastolic_myo_shapes(n,:) = all_training_diastolic_myo_shapes(c,:);
        %             DETERMINE_systolic_myo_shapes(n,:) = all_training_systolic_myo_shapes(c,:);
        %
        %             DETERMINE_dia_sys_myo_shapes(n,:) = all_training_dia_sys_myo_shapes(c,:);
        %
        %         end
        
    end
    
    dia_shapes = all_training_diastolic_myo_shapes;
    sys_shapes = all_training_diastolic_myo_shapes;
    dia_sys_shapes = all_training_dia_sys_myo_shapes;
    
    dia_shapes( ~any(dia_shapes,2), : ) = [];
    sys_shapes( ~any(sys_shapes,2), : ) = [];
    dia_sys_shapes( ~any(dia_sys_shapes,2), : ) = [];
    
    dia_myo_mean = mean(dia_shapes);
    sys_myo_mean = mean(sys_shapes );
    dia_sys_myo_mean = mean(dia_sys_shapes) ;
    
    %calculate means
    %     dia_myo_mean = diastolic_myo_shapes(1).xyz;
    %     sys_myo_mean = systolic_myo_shapes(1).xyz;
    %     for iS = 2:size(diastolic_myo_shapes,2)
    %         dia_myo_mean = dia_myo_mean + diastolic_myo_shapes(iS).xyz;
    %         sys_myo_mean = sys_myo_mean + systolic_myo_shapes(iS).xyz;
    %     end
    %     dia_myo_mean = dia_myo_mean./(size(diastolic_myo_shapes,2)-1);
    %     sys_myo_mean = sys_myo_mean./(size(systolic_myo_shapes,2)-1);
    
    % % update the reference shapes for the next round of procrustes
    % dia_endo_reference = dia_endo_mean;
    % dia_epi_reference = dia_epi_mean;
    % sys_endo_reference = sys_endo_mean;
    % sys_epi_reference = sys_epi_mean;
    
    dia_myo_reference = reshape(dia_myo_mean,[2*1089 3]);
    sys_myo_reference = reshape(sys_myo_mean,[2*1089 3]);
    dia_sys_reference = reshape(dia_sys_myo_mean, [4*1089 3]);
    
    procrustes_iteration = p
end



disp('finished procrustes')
end