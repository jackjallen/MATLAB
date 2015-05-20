%% VG code for visualising the project data
% recieved: 20/04/15
close all
%%
%VG
%pat1_endo=xlsread('F:\Images\STACOM\stacom-training-pack-clean\output-stacom-clean\SSM0001.ED.endo.vertices.csv');
pat1_endo=xlsread('C:\Users\jesu2687\Documents\MATLAB\ErnestoCode\Subject1\SSM0001.ES.endo.vertices.csv');
figure
plot3(pat1_endo(:,1),pat1_endo(:,2),pat1_endo(:,3),'o')
axis equal; axis tight;
%VG
%pat1_epi=xlsread('F:\Images\STACOM\stacom-training-pack-clean\output-stacom-clean\SSM0001.ED.epi.vertices.csv');
pat1_epi=xlsread('C:\Users\jesu2687\Documents\MATLAB\ErnestoCode\Subject1\SSM0001.ES.epi.vertices.csv');
hold on
plot3(pat1_epi(:,1),pat1_epi(:,2),pat1_epi(:,3),'or')
title 'SSM0001 endocardium'
legend 'diastolic' 'systolic'
%VG
%trifac=xlsread('F:\Images\STACOM\stacom-training-pack-clean\TriangleFaces.csv');
trifac=xlsread('C:\Users\jesu2687\Documents\MATLAB\ErnestoCode\Subject1\TriangleFaces.csv');
p.vertices=pat1_endo;
p.faces=trifac;
figure
axis equal; axis tight;
patch(p,'EdgeColor','g','FaceColor','red');
lighting gouraud
camlight headlight
hold on
patch('vertices', pat1_epi, 'faces', p.faces,'edgecolor','k','facecolor','g')
patch('vertices', pat1_endo, 'faces', p.faces,'edgecolor','k','facecolor','b')
%% JA
% 20/04/15 triangle areas
% area = 0.5 * base *height

%plot3(trifac(:,1),trifac(:,2),trifac(:,3))
%plot3(pat1_endo(:,1),pat1_endo(:,2),pat1_endo(:,3))



%%
%VG
%pat1_endo=xlsread('F:\Images\STACOM\stacom-training-pack-clean\output-stacom-clean\SSM0001.ED.endo.vertices.csv');
pat1_endo=xlsread('C:\Users\jesu2687\Documents\stacom-training-pack-clean\stacom-training-pack-clean\output-stacom-clean\output-stacom-clean\SSM0002.ED.endo.vertices.csv');
figure
plot3(pat1_endo(:,1),pat1_endo(:,2),pat1_endo(:,3),'o')
axis equal; axis tight;
%VG
%pat1_epi=xlsread('F:\Images\STACOM\stacom-training-pack-clean\output-stacom-clean\SSM0001.ED.epi.vertices.csv');
pat1_epi=xlsread('C:\Users\jesu2687\Documents\stacom-training-pack-clean\stacom-training-pack-clean\output-stacom-clean\output-stacom-clean\SSM0002.ED.epi.vertices.csv');
hold on
plot3(pat1_epi(:,1),pat1_epi(:,2),pat1_epi(:,3),'or')
%VG
%trifac=xlsread('F:\Images\STACOM\stacom-training-pack-clean\TriangleFaces.csv');
trifac=xlsread('C:\Users\jesu2687\Documents\stacom-training-pack-clean\stacom-training-pack-clean\TriangleFaces.csv');
p.vertices=pat1_endo;
p.faces=trifac;
figure
axis equal; axis tight;
patch(p,'EdgeColor','g','FaceColor','red');
lighting gouraud
camlight headlight