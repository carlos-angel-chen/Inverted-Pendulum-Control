%% State feedback

clear all; close all; clc

%% Constantes del modelo.
m = 0.1818;
M = 0.6042;
Mt = 10;
l = 0.09; %cm 
g = 10;
J = 0.45e-3;
%JJt = J+m*l^2;
Jt = 2.744e-3;
mu = Mt*Jt-(m^2)*(l^2);

%% Modelo sistema : Matrices (A,B)


%% Realimentación de estados
%x = [p,q,v,w]
K = 1.0e+03 * [-8.4073,3.7519,-2.2419,0.4756]
