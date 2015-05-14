%% myocardium visualisation
close all
%visualise the endo and epi edge points (labelled), with all the other points from
%endo and epi.
figure
%subplot 221
hold on
plot3(data(i).diastolic.myo.B.xyz(:,1),data(i).diastolic.myo.B.xyz(:,2), data(i).diastolic.myo.B.xyz(:,3))
text(data(i).diastolic.myo.B.xyz(:,1) ,  data(i).diastolic.myo.B.xyz(:,2) ,  data(i).diastolic.myo.B.xyz(:,3) , arrayfun( @(id)sprintf('%d',id) , 1:size(data(i).diastolic.myo.B.xyz,1) , 'un',false ) )
plot3(data(i).diastolic.full.xyz(:,1),data(i).diastolic.full.xyz(:,2), data(i).diastolic.full.xyz(:,3),'o')
legend('endo and epi edge points (myo B)', 'all endo and epi points')

%visualise surfaces: endo, epi and myo lid
%cla
figure
%subplot 222
patch('vertices',data(i).diastolic.endo.xyz,'faces',data(i).diastolic.endo.tri,'facecolor','none','EdgeColor','red')
patch('vertices',data(i).diastolic.epi.xyz,'faces',data(i).diastolic.epi.tri,'facecolor','none','EdgeColor','blue')
patch('vertices',data(i).diastolic.myo.B.xyz,'faces',myoB.tri,'facecolor','green')
legend('endo', 'epi', 'myo lid')

% visualise the myocardium mesh (solid)
%cla
figure
%subplot 223
patch('vertices',data(i).diastolic.myo.xyz,'faces',data(i).diastolic.myo.tri,'facecolor','red')
legend('full myocardium')

figure
%subplot 223
patch('vertices',data(i).systolic.myo.xyz,'faces',data(i).systolic.myo.tri,'facecolor','red')
legend('full myocardium')

%visualise center of mass within myocardium mesh
%cla
figure
%subplot 224
patch('vertices',data(i).diastolic.myo.xyz,'faces',data(i).diastolic.myo.tri,'facecolor','g','facealpha',0.1);
hold on;
plot3( data(i).diastolic.myoCenterOfMass(1) ,  data(i).diastolic.myoCenterOfMass(2) ,  data(i).diastolic.myoCenterOfMass(3) , '*r','markers',20 ); hold off
legend('full myocardium', 'myocardium center of mass')
