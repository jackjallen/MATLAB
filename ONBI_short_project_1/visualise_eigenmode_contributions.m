%% visualise eigenvalue contributions
dia_myo_normalisedEigvalues = sorted_dia_myo_eigenvalues/(sum(sorted_dia_myo_eigenvalues));
sys_myo_normalisedEigvalues = sorted_sys_myo_eigenvalues/(sum(sorted_sys_myo_eigenvalues));
dia_sys_myo_normalisedEigvalues = sorted_dia_sys_myo_eigenvalues/(sum(sorted_dia_sys_myo_eigenvalues));
dia_sys_CS = cumsum(dia_sys_myo_normalisedEigvalues);
dia_CS = cumsum(dia_myo_normalisedEigvalues);
sys_CS = cumsum(sys_myo_normalisedEigvalues);

figure
nEig = 15;
hold on
plot(dia_CS(1:100))
plot(sys_CS(1:100))
for eigVal = 1:nEig
    plot(eigVal,(dia_CS(eigVal)),'o')
    dia_CS(eigVal)
end
for eigVal = 1:nEig
    plot(eigVal,(sys_CS(eigVal)),'*')
end
xlabel 'eigenmode'
ylabel 'contribution'
legend 'diastole' 'systole'
hold off


figure
hold on
plot(dia_sys_CS(1:100))
for eigVal = 1:nEig
    plot(eigVal,(dia_sys_CS(eigVal)),'o')
    dia_sys_CS(eigVal)
end
xlabel 'eigenmode'
ylabel 'contribution'
legend 'diastole and systole'
hold off

