function [area] = calcTriMeshArea(sides)
% Calculation of the area of a triangle using heron's forumla.
% Receives a matrix of side lengths.
% Returns total area as a scalar

% semiperimetre s
s = (sides(:,1) + sides(:,2) + sides(:,3))./2;
% total area (scalar)
area = sum(sqrt(s.*(s-sides(:,1)).*(s-sides(:,2)).*(s-sides(:,3))));
end