function [ x ] = flsa_SLEP( y, lambda1, lambda2 )
%FLSA Summary of this function goes here
%   Detailed explanation goes here

n = length(y);
nn = n - 1;

% the duality gap for termination
tol = 1e-10;

% the maximal number of iterations
maxStep = 100;

% the starting point
z0 = zeros(nn, 1);

[x, z, infor] = flsa(y, z0, lambda1, lambda2, n, maxStep, tol, 1, 6);

end