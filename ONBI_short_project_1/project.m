%% ONBI Short project 1
% Jack Allen
% Supervisor: Vicente Grau
%
clear all
close all
clc
%%
% read data (adapted from VG code)
%patient 1, endocardium
pat1_endo=xlsread('C:\Users\jesu2687\Documents\stacom-training-pack-clean\stacom-training-pack-clean\output-stacom-clean\output-stacom-clean\SSM0001.ED.endo.vertices.csv');
endo = pat1_endo;
%epicardium
pat1_epi=xlsread('C:\Users\jesu2687\Documents\stacom-training-pack-clean\stacom-training-pack-clean\output-stacom-clean\output-stacom-clean\SSM0001.ED.epi.vertices.csv');
epi = pat1_epi;
%triangle faces
trifac=xlsread('C:\Users\jesu2687\Documents\stacom-training-pack-clean\stacom-training-pack-clean\TriangleFaces.csv');
p.vertices=endo;
p.faces=trifac;

%% calculate mean shape

%load the data points
load('project1_data.mat');
ED_endo = cell(1,400);
tic;
for i = 1:9
    ED_endo{i} = xlsread(['C:\Users\jesu2687\Documents\stacom-training-pack-clean\stacom-training-pack-clean\output-stacom-clean\output-stacom-clean\SSM000' num2str(i) '.ED.endo.vertices.csv']);
    ES_endo{i} = xlsread(['C:\Users\jesu2687\Documents\stacom-training-pack-clean\stacom-training-pack-clean\output-stacom-clean\output-stacom-clean\SSM000' num2str(i) '.ES.endo.vertices.csv']);
    ED_epi{i} = xlsread(['C:\Users\jesu2687\Documents\stacom-training-pack-clean\stacom-training-pack-clean\output-stacom-clean\output-stacom-clean\SSM000' num2str(i) '.ED.epi.vertices.csv']);
    ES_epi{i} = xlsread(['C:\Users\jesu2687\Documents\stacom-training-pack-clean\stacom-training-pack-clean\output-stacom-clean\output-stacom-clean\SSM000' num2str(i) '.ES.epi.vertices.csv']);
    
    i    
end
for i = 10:99
    ED_endo{i} = xlsread(['C:\Users\jesu2687\Documents\stacom-training-pack-clean\stacom-training-pack-clean\output-stacom-clean\output-stacom-clean\SSM00' num2str(i) '.ED.endo.vertices.csv']);
    ES_endo{i} = xlsread(['C:\Users\jesu2687\Documents\stacom-training-pack-clean\stacom-training-pack-clean\output-stacom-clean\output-stacom-clean\SSM00' num2str(i) '.ES.endo.vertices.csv']);
    ED_epi{i} = xlsread(['C:\Users\jesu2687\Documents\stacom-training-pack-clean\stacom-training-pack-clean\output-stacom-clean\output-stacom-clean\SSM00' num2str(i) '.ED.epi.vertices.csv']);
    ES_epi{i} = xlsread(['C:\Users\jesu2687\Documents\stacom-training-pack-clean\stacom-training-pack-clean\output-stacom-clean\output-stacom-clean\SSM00' num2str(i) '.ES.epi.vertices.csv']);
    
    i    
end
for i = 100:400
    ED_endo{i} = xlsread(['C:\Users\jesu2687\Documents\stacom-training-pack-clean\stacom-training-pack-clean\output-stacom-clean\output-stacom-clean\SSM0' num2str(i) '.ED.endo.vertices.csv']);
    ES_endo{i} = xlsread(['C:\Users\jesu2687\Documents\stacom-training-pack-clean\stacom-training-pack-clean\output-stacom-clean\output-stacom-clean\SSM0' num2str(i) '.ES.endo.vertices.csv']);
    ED_epi{i} = xlsread(['C:\Users\jesu2687\Documents\stacom-training-pack-clean\stacom-training-pack-clean\output-stacom-clean\output-stacom-clean\SSM0' num2str(i) '.ED.epi.vertices.csv']);
    ES_epi{i} = xlsread(['C:\Users\jesu2687\Documents\stacom-training-pack-clean\stacom-training-pack-clean\output-stacom-clean\output-stacom-clean\SSM0' num2str(i) '.ES.epi.vertices.csv']);
    
    i
end
toc;
%findmean
ED_endo_sum = ED_endo{1};
tic;
for i = 2:400
tmp = ED_endo{1,1};
ED_endo_sum = ED_endo_sum + tmp;
end
toc;
ED_endo_mean = ED_endo_sum/400;
%% Covariance matrix

%% calculate triangle side lengths
endo_sides = calcTriSides(trifac, endo);
epi_sides = calcTriSides(trifac, endo);
%% calculate total areas(heron's forumla)
endo_area = calcTriMeshArea(endo_sides);
epi_area = calcTriMeshArea(epi_sides);
%% calculate coordinates of the centroids
endo_centroid = mean(endo);
epi_centroid = mean(epi);
