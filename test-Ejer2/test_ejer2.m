%% State feedback

clear all; close all; clc

%% Constantes del modelo.
m = 0.1818; 
M = 0.6042;
Mt = m+M;
l = 0.09; %cm 
g = 9.81;
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

% Funcion de transferencia del modelo linealizado
T = zpk(Gss)
stepinfo(T);

%% C2_1 - Loop shaping: INTERNO con T2
s = tf('s');
T2 = 1/(40*10*2*pi);   %Arranco considerando un periodo de muestreo chico, parto de una situación mas ideal. 
Pade = (1-s*T2/4)/(1+s*T2/4);
C2_1 = 1;
Pap2 = 8.6614*(s+8.172)/(s-8.172);
Pmp2 = 1/((s+8.172)^2);
L2_1 = Pap2*Pmp2*C2_1*Pade
figure('Name', 'Bode2_1'); bode(L2_1)
margin(L2_1)
grid on

%% C2_2 - Proponemos un control PD con el objetivo de estabilizar la salida
% Proponemos una K y un alpha que multiplica al control para obtener un pm 
% igual a 60
K2 = 1;
alpha = 1;
C2_2 = K2*(s+8.172)/(s+alpha*8.172);
L2_2 = Pap2*Pmp2*Pade*C2_2
figure('Name', 'Bode2_1'); bode(L2_2)
margin(L2_2)
grid on
T2_2 = L2_2/(1+L2_2);
figure('Name', 'Step2_2'); step(T2_2);
grid on;
CS2_2 = C2_2/(1+L2_2)
figure('Name', 'StepCS2_3'); step(CS2_2)

%% C2_3 - Configuro los valores de K y alpha para tener una salida estable 
% y pm = 60
K2 = db2mag(55);
alpha = 15;
C2_3 = K2*(s+8.172)/(s+alpha*8.172);
L2_3 = Pap2*Pmp2*Pade*C2_3
figure('Name', 'Bode2_1'); bode(L2_3)
margin(L2_3);
grid on
T2_3 = L2_3/(1+L2_3);
figure('Name', 'Step2_2'); step(T2_3);
grid on;
CS2_3 = C2_3/(1+L2_3)
figure('Name', 'StepCS2_3'); step(CS2_3)

%% Hallar la planta que ve C1 (Reurison)
s  = tf('s');
C2 = C2_3; %controlador C2 (interno)
C2    = ss(C2);
C2.u  = 'tita';
C2.y  = 'u2';
Gss   = ss(A,B,C,0); %Sistema completo
Gss.u = 'u'; %soma de u(C1) + u(C2)
Gss.y = 'y';
Sum = sumblk('u = u1 - u2');

Sysq   = ss([],[],[],[1 0]); %q
Sysq.u = 'y';
Sysq.y = 'q';

Systita   = ss([],[],[],[0 1]); %tita
Systita.u = 'y';
Systita.y = 'tita';

% Transferencia "vista" por C1.
GC1 = connect(Gss, C2, Sum, Sysq, Systita, 'u1', 'q');
GC1 = zpk(GC1);

%% C1_1 - Analizando la planta que ve el C1
C1_1 = 1;
P1 = GC1;
pap = (7.648-s)/(s+7.648);
pmp = -(1.4526*(s+122.6)*(s+7.648)^2)/(s^2*(s+8.172)*(s^2+114.4*s+3869));

L1_1 = C1_1*pap*pmp*Pade;
figure('Name', 'Bode1_1'); bode(L1_1)
margin(L1_1)
grid on

T1_1 = L1_1/(1+L1_1);
figure('Name', "Step"), step(T1_1);
grid on;

%% C1_2 - Agrego accion integral y elimino polos estables
% Controlador impropio
K1 = db2mag(6.1);
C1_2 = (1/s)*(1/pmp);

L1_2 = K1*C1_2*pap*pmp*Pade;
figure('Name', 'Bode1_2'); bode(L1_2)
margin(L1_2)
grid on

T1_2 = L1_2/(1+L1_2);
figure('Name', "Step"), step(T1_2);
grid on;

%% C1_3 - Agrego pendientes en altas y bajas frecuencias
K1 = db2mag(53.9);
C1_3 = (1/s)*(1/pmp)*(1/(s+20)^2)*((s+0.1)/(s));

L1_3 = K1*C1_3*pap*pmp*Pade;
figure('Name', 'Bode1_3'); bode(L1_3)
margin(L1_3)
grid on

T1_3 = L1_3/(1+L1_3);
figure('Name', "Step"), step(T1_3);
grid on;

CS1_3 = C1_3/(1+L1_3);
figure('Name', 'StepCS2_3'); step(CS1_3);
grid on;

% Error en estado estacionario
1 - dcgain(T1_3)

%% Graficos
%Lazo Abierto
L = L1_3;
subplot(3,2,1)
bode(L)
title("Bode de L(s) lazo abierto")
margin(L)
grid on

%Lazo Cerrado
S  = 1/(1+L); 
T  = L/(1+L);
CS = C1_3/(1+L);
PS = P1/(1+L);

%Graficos bode
subplot(3,2,2)
bode(S,T,CS,PS)
title("Bode de S(s),T(s),CS(s) y PS(s)")
grid on

%Accion de control
subplot(3,2,3)
step(CS)
title("Accion de Control C(s)S(s)")
grid on
AccionDeControl_info = stepinfo(CS)

%Respuesta temporal
subplot(3,2,4)
step(PS)
title("Respuesta Temporal P(s)S(s)")
grid on
RespuestaPS_info = stepinfo(T)

%Respuesta temporal
subplot(3,2,5)
step(T)
title("Respuesta Temporal T(s)")
grid on
RespuestaTemporal_info = stepinfo(T)




