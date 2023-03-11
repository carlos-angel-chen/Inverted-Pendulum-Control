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

%Variables de estado y entrada
syms x1 x2 x3 x4;       %q theta q_d theta_d
syms u;

%Ecuaciones no lineales. Modelo de Astrom & Murray 3-13.

x1d = x3;
x2d = x4;
x3d = (-m*l*sin(x2)*(x4^2)+m*g*(m*(l^2)/Jt)*sin(x2)*cos(x2)+u)/(Mt-m*(m*(l^2)/Jt)*(cos(x2))^2);
x4d = (-m*(l^2)*sin(x2)*cos(x2)*(x4^2)+Mt*g*l*sin(x2)+l*cos(x2)*u)/((Jt*Mt/m)-m*(l*cos(x2))^2);

h1 = [x1; x2];

% Vectores
x = [x1; x2; x3; x4];
f = [x1d; x2d; x3d; x4d];
h = [h1];

%Linealizacion del sistema (Jacobiano)
A = jacobian(f, x);  %df/dx
B = jacobian(f, u);  %df/du
C = jacobian(h, x);  %dh/dx
D = jacobian(h, u);  %dh/du

% Puntos de equlibrio
x1_eq = 0;
x2_eq = 0;
x3_eq = 0;
x4_eq = 0;
u_eq = 0;

Xeq = [x1_eq x2_eq x3_eq x4_eq];
delta_xeq = [0 0 0 0];
Ueq = [u_eq];

%Calculo mmatrices en pto de equlibrio
A = subs(A, {'x1','x2','x3', 'x4', 'u'}, {x1_eq, x2_eq, x3_eq, x4_eq, u_eq});
B = subs(B, {'x1','x2','x3', 'x4','u'}, {x1_eq, x2_eq, x3_eq, x4_eq, u_eq});
C = subs(C, {'x1','x2','x3', 'x4','u'}, {x1_eq, x2_eq, x3_eq, x4_eq, u_eq});
D = subs(D, {'x1','x2','x3', 'x4','u'}, {x1_eq, x2_eq, x3_eq, x4_eq, u_eq});

% Gss : Sistema linealizado en espacio de estados
A = double(A);
B = double(B);
C = double(C);
D = double(D);

%% Estabilidad del sistema linealizado
% Analizar la estabilidad del sistema linealizado
Gss = ss(A,B,C,D);
eig(A)

%% Funcion de transferencia del modelo linealizado
T = zpk(Gss);
stepinfo(T);

%% Planta del angulo
s = tf('s');
P_ang = 0.60216/((s-7.76)*(s+7.76));
Ts = 1/(5*10*2*pi);
pade = (1 - s*Ts/4)/(1 + s*Ts/4);
pap = pade*(s+7.76)/(s-7.76);
pmp = 1/((s+7.76)^2);
P = pap*pmp;















