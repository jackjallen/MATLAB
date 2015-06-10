%% visualise training data

bestParam3 = 1 %EF
b = [bestParam1(1);bestParam2(1);bestParam3]';
% b = [35;36;37];
figure
hold on
%histogram(trainingdata(1:100,b))
plot3D(trainingdata(1:100,b),'o')
plot3D(trainingdata(101:200,b),'o')
hold off
legend ' DETERMINE' 'MESA'
xlabel ([ num2str(bestParam1(1)) ])
ylabel ([ num2str(bestParam2(1)) ])
zlabel ([ num2str(bestParam3) ' (EF) '])