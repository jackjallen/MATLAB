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

%allocate matrix to store side lengths of triangles

%% calculate triangle side lengths
endo_sides = calcTriSides(faces, vertices);
epi_sides = calcTriSides(faces, vertices);
%% Area calculation (heron's forumla)
% semiperimetre s
s = (sides(:,1) + sides(:,2) + sides(:,3))./2;
% list of triangle areas
areas = sqrt(s.*(s-sides(:,1)).*(s-sides(:,2)).*(s-sides(:,3)));
total_area = sum(areas);

%% Centroids
%coordinates of centroids
endo_centroid = mean(endo);
epi_centroid = mean(epi);
