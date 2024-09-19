function [GEOMETRY, Geo_id] = mesh2GEOMETRY(GCOORD, ELEM2NODE, Point_id)
% map triangulated mesh (by external mesh generator triangle) into an
% updated GEOMETRY containing all new triangle vertices added by the
% mesh generator
%
% This function works only if input to triangle is hierarchical, with
% higher order interfaces having a preference over lower order (as done by 
% generate_meshGEO)
%
% Author: Javier García-Pintado, 2020-03

outputArg1 = inputArg1;
outputArg2 = inputArg2;
end

