%% connect 

s  = tf('s');
C2 = ; %controlador C2 (interno)
C2    = ss(C2);
C2.u  = 'tita';
C2.y  = 'u2';
Gss   = ss(A,B,C,0); %Sistema completo
Gss.u = 'u'; %soma de u(C1) + u(C2)
Gss.y = 'y';
Sum = sumblk('u = u1 + u2');

Sysq   = ss([],[],[],[1 0]); %q
Sysq.u = 'y';
Sysq.y = 'q';

Systita   = ss([],[],[],[0 1]); %tita
Systita.u = 'y';
Systita.y = 'tita';

% Transferencia "vista" por C1.
GC1 = connect(Gss, C2, Sum, Sysq, Systita, 'u1', 'q');
GC1 = zpk(GC1);