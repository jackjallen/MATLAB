%% FIND MYOCARDIUM THICKNESS EIGENMODES

for i = 1:401
dia_sys_dEPI2ENDOs(i,:) = [dia_dEPI2ENDOs(i,:)' ; sys_dEPI2ENDOs(i,:)']' ;
for m = data(1).MESA_indices
    MESA_dia_sys_dEPI2ENDOs(m,:) = [dia_dEPI2ENDOs(m,:)' ; sys_dEPI2ENDOs(m,:)']' ;
end
end
%%  PCA on myocardium thicknesses

dia_sys_dEPI2ENDOs_cov_mat = cov(dia_sys_dEPI2ENDOs);

%% PCA - find eigenvectors of covariance matrix 
[dia_sys_dEPI2ENDOs_eigenvectors,dia_sys_dEPI2ENDOs_eigenvalues]=eig(dia_sys_dEPI2ENDOs_cov_mat);

nEigenvalues = 15;
sorted_dia_sys_dEPI2ENDOs_eigenvalues = sort(dia_sys_dEPI2ENDOs_eigenvalues(:), 'descend');

principle_dia_sys_dEPI2ENDOs_eigenvalues = sorted_dia_sys_dEPI2ENDOs_eigenvalues(1:nEigenvalues);

sorted_dia_sys_dEPI2ENDOs_eigenvalues_contributions = 100*(sorted_dia_sys_dEPI2ENDOs_eigenvalues./(sum(sorted_dia_sys_dEPI2ENDOs_eigenvalues)));

%% PCA - contribution from the chosen eigenvectors, as a percentage
principle_dia_sys_dEPI2ENDOs_eigenvalues_contribution = sorted_dia_sys_dEPI2ENDOs_eigenvalues_contributions(1:nEigenvalues);

principle_dia_sys_dEPI2ENDOs_eigenvectors = zeros(2178, nEigenvalues);

for n = 1:nEigenvalues  
    [dia_sys_dEPI2ENDOs_eRows(1,n), dia_sys_dEPI2ENDOs_eCols(1,n)] = find(dia_sys_dEPI2ENDOs_eigenvalues == principle_dia_sys_dEPI2ENDOs_eigenvalues(n));
  
    principle_dia_sys_dEPI2ENDOs_eigenvectors(:,n) = dia_sys_dEPI2ENDOs_eigenvectors(:, dia_sys_dEPI2ENDOs_eCols(1,n)); 
end

%% Find max model parameter value (b) , using +/- sqrt(eigenvalue) as b
% b can be increased later.

dia_sys_dEPI2ENDOs_max_b = sqrt(principle_dia_sys_dEPI2ENDOs_eigenvalues);


%% PDM - find mean shape
dia_sys_dEPI2ENDOs_mean = mean(dia_sys_dEPI2ENDOs);

%% PDM - find b values for the DETERMINE and MESA cases
for b = 1:nEigenvalues
    for i = data(1).DETERMINE_indices'
        for d = find(data(1).DETERMINE_indices==i)
            
           DETERMINE_dia_sys_dEPI2ENDOs_b(d,b) = (principle_dia_sys_dEPI2ENDOs_eigenvectors(:,b)')*(dia_sys_dEPI2ENDOs(i,:)' - dia_sys_dEPI2ENDOs_mean');
            
        end
    end
    
    for i = data(1).MESA_indices'
        for d = find(data(1).MESA_indices==i)
            %         d
            MESA_dia_sys_dEPI2ENDOs_b(d,b) = (principle_dia_sys_dEPI2ENDOs_eigenvectors(:,b)')*(dia_sys_dEPI2ENDOs(i,:)' - dia_sys_dEPI2ENDOs_mean');
            
        end
    end
end


%% visualise eigenvalue contributions
dia_sys_dEPI2ENDOs_normalisedEigvalues = sorted_dia_sys_dEPI2ENDOs_eigenvalues/(sum(sorted_dia_sys_dEPI2ENDOs_eigenvalues));


dia_sys_CS = cumsum(dia_sys_dEPI2ENDOs_normalisedEigvalues);

figure
nEig = 15;
hold on
% plot(dia_CS(1:100))
plot(dia_sys_CS(1:100))
% for eigVal = 1:nEig
%     plot(eigVal,(dia_CS(eigVal)),'o')
%     dia_CS(eigVal)
% end
for eigVal = 1:nEig
    plot(eigVal,(dia_sys_CS(eigVal)),'*')
end
xlabel 'eigenmode'
ylabel 'contribution'
legend 'diastole' 'systole'
hold off

%% visualise modes of variation - static
% close all

x = [-60 ; 40];
y = [-65 ;50];
z = [-110 ; 20];
% figure('name','diastolic modes of variation')
run('visualisemodes_diastole_myo_thicknesses_bullseyes.m')
% run('visualisemodes_diastole__combined_modes_myo_thicknesses.m')
% figure('name','systolic modes of variation')
run('visualisemodes_systole_myo_thicknesses_bullseyes.m')
% run('visualisemodes_systole__combined_modes_myo_thicknesses.m')
