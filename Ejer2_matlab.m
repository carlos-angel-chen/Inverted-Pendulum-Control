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
T1 = 10   %Arranco considerando un periodo de muestreo chico, parto de una situación mas ideal. 
Pade = (1-s*T1/4)/(1+s*T1/4);
Pap_i = Pade*(0.60216*(s+7.76)/(s-7.76))
Pmp_i = 1/(s+7.76)^2
C1_i = 1
L1_i = Pap_i*Pmp_i*C1_i
figure('Name', 'Bode1'); bode(L1_i)
Tf1_i = L1/(1+L1_i)
figure('Name', 'Step1'); step(Tf1_i)
CS1_i = C1/(1+L1_i)
figure('Name', 'StepCS1'); step(CS1_i)

%% Agrego un compensador Lead para agregar mas fase a bajas frecuencias
%Agrego un polo en alta frecuencia para eliminar ruido. Aprox 20 veces ubicacion de frecuencia de cruce (1.86Hz). 
%Agrego un lag antes de frecuencia de cruce para tener mas pendiente.
%Modifico la ganancia y voy jugando con el Ts hasta que llego a un sistema que cumpla especificaciones. 

K2_i = db2mag(53) %Agrego ganancia para MF = 60°
C2_i = K2_i*((s+7.76)^2*(s+0.8/10))/(s^2 * (s+0.8/100) * (s+1*20)^2);
L2_i = Pap*Pmp*C2
figure('Name', 'Bode2'); bode(L2_i)
Tf2_i = L2/(1+L2_i)
figure('Name', 'Step2'); step(Tf2_i)
CS2_i = C2/(1+L2_i)
figure('Name', 'StepCS2'); step(CS2_i)

%% Loop shaping: EXTERNO con T1
s = tf('s');
Pap_e = Pade*(0.10099*(s+7.76)*(s-7.722))/((s-7.76)*(s+7.722))
Pmp_e = (s+7.722)/(s+7.76)^2
C1_e = 1
L1_e = Pap_e*Pmp_e*C1_e
figure('Name', 'Bode1'); bode(L1_e)
Tf1_e = L1_e/(1+L1_e)
figure('Name', 'Step1'); step(Tf1_e)
CS1_e = C1_e/(1+L1_e)
figure('Name', 'StepCS1'); step(CS1_e)

%% Agrego un compensador Lead para agregar mas fase a bajas frecuencias FALTA TERMINAR DE CORREGIR CONTROLADOR
%Agrego un polo en alta frecuencia para eliminar ruido. Aprox 20 veces ubicacion de frecuencia de cruce (1.86Hz). 
%Agrego un lag antes de frecuencia de cruce para tener mas pendiente.
%Modifico la ganancia y voy jugando con el Ts hasta que llego a un sistema que cumpla especificaciones. 

K2_e = db2mag(0.722) %Agrego ganancia para MF = 60°
C2_e = K2_e*((s+7.76))/(s);
L2_e = Pap_e*Pmp_e*C2_e
figure('Name', 'Bode2'); bode(L2_e)
Tf2_e = L2_e/(1+L2_e)
figure('Name', 'Step2'); step(Tf2_e)
CS2_e = C2_e/(1+L2_e)
figure('Name', 'StepCS2'); step(CS2_e)




