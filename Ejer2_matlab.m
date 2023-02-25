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
T = zpk(Gss)
stepinfo(T);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Loop shaping: INTERNO con T2
s = tf('s');
T2 = 1/(40*10*2*pi);   %Arranco considerando un periodo de muestreo chico, parto de una situación mas ideal. 
Pade = (1-s*T2/4)/(1+s*T2/4);
C2_1 = 1;
Pap2 = 0.60216*(s+7.76)/(s-7.76);
Pmp2 = 1/((s+7.76)^2);
L2_1 = Pap2*Pmp2*C2_1*Pade
figure('Name', 'Bode2_1'); bode(L2_1)
margin(L2_1)
grid on

%% 
C2_2 = db2mag(34.7)*(1/Pmp2)*(1/s);
L2_2 = Pap2*Pmp2*Pade*C2_2
figure('Name', 'Bode2_1'); bode(L2_2)
margin(L2_2)
grid on

%% 
C2_3 = db2mag(37.5)*(1/Pmp2)*(1/s)*(1/(s/(30*20)+1)^2)*((s+1)/(s+0.1));
L2_3 = Pap2*Pmp2*Pade*C2_3
figure('Name', 'Bode2_1'); bode(L2_3)
margin(L2_3)
grid on

Tf2_3 = L2_3/(1+L2_3)
figure('Name', 'Step2_3'); step(Tf2_3)
CS2_3 = C2_3/(1+L2_3)
figure('Name', 'StepCS2_3'); step(CS2_3)

%%
nyqlog(L2_3)

%% Loop shaping: COMPLETO
C1_1 = ;
L1_1 = Pap2*Pmp2*Pade*C1_1
figure('Name', 'Bode1_1'); bode(L1_1)
margin(L1_1)
grid on

Tf1_1 = L1_1/(1+L1_1)
figure('Name', 'Step1_1'); step(Tf1_1)
CS1_1 = C1_1/(1+L1_1)
figure('Name', 'StepCS1_1'); step(CS1_1)



