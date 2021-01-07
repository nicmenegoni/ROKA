function [d] = dist2points(P1,P2)
%% INFORMATION
% This function calculates the distance between two 3D points. This script
% was tested to respect the 'norm' function of Matlab and it is around 16%
% more faster.
% Input parametrs:
%   P1 = vector containing the x, y and z coordinates of the point 1
%   P2 = vector containing the x, y and z coordinates of the point 2 
% Output parameters:
%   d = sistance between point 1 and point 2
%% SCRIPT
d = sqrt( (P1(1)-P2(1))^2 + (P1(2)-P2(2))^2 + (P1(3)-P2(3))^2);
end