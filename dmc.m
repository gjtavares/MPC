close all; clc; clear all

%% Funcao de transferencia em Laplace

NumGs = 3.2132e-04;
DenGs = [3.4051 3.2132e-04];
Gs = tf(NumGs, DenGs);
stepinfo(Gs) % ts = 4.1457e+04

NumP1s = -0.0221;
DenP1s = [3.4051 3.2132e-04];
P1s = tf(NumP1s, DenP1s);

NumP2s = 2.3923e-07;
DenP2s = [3.4051 3.2132e-04];
P2s = tf(NumP2s, DenP2s);
%% Período de amostragem

Ts = 1382; % ts/15 <= Ts <= ts/6,  Ts = (ts/2)/15, ts/2 --> Tempo de acomodação em malha fechada

%% Discretizacao do modelo em tempo continuo:

Gz = c2d(Gs, Ts, 'zoh')
P1z = c2d(P1s, Ts, 'zoh')
P2z = c2d(P2s, Ts, 'zoh')

%Parametros

Ph = 10; % Horizonte de predição
Mh = 5; % Horizonte de controle
salto = 30;
lambda = 1;

%% Resposta ao Degrau Gz

NumGz = [0 0.1223];
DenGz = [1 -0.8777];
atraso = 0;

gd = resposta_degrau(NumGz,DenGz,atraso,salto);
gd((31):(40),1) = gd(30,1);

%% Resposta ao degrau Disturbio

NumP1z = [0 -8.409];
DenP1z = [1 -0.8777];

dg = resposta_degrau(NumP1z,DenP1z,atraso,salto);
dg(31:40,1) = dg(30,1);

%% Matriz dinâmica

G = zeros(Ph,Mh);

G(1:10,1) = gd(1:10,1);
G(2:10,2) = G(1:9,1);
G(3:10,3) = G(1:8,1);
G(4:10,4) = G(1:7,1);
G(5:10,5) = G(1:6,1);

%% Matriz D

D = zeros(Ph,Mh);

D(1:10,1) = dg(1:10,1);
D(2:10,2) = D(1:9,1);
D(3:10,3) = D(1:8,1);
D(4:10,4) = D(1:7,1);
D(5:10,5) = D(1:6,1);

%% Sequencia de Controle

R=chol(G'*G+lambda*eye(size(G'*G)));
K=R\(R'\G');

%% 

per_amost = Ts;
t_final = Ts*120;
N = t_final/per_amost; %N° amostras
tempo = 0:per_amost:t_final-per_amost; %Vetor de tempo

w=35*ones(N+Ph,1);
wp=w(1:N,1);
p=0.5*[zeros(30,1);ones(50,1);zeros(100,1)];
pp=p(1:N,1);

%% Definação das variáveis

delta_u(salto,1)=0;
delta_p(salto,1)=0;
deltap1=0;

y_saida(N,1)=0;
u_atual=0;

Denconv = conv(DenGz,DenP1z);
Uconv = conv(NumGz,DenP1z);
Pconv = conv(NumP1z,DenGz);

na = size(Denconv,2);
nb = size(Uconv,2);
nbp = size(Pconv,2);

yna = zeros(na,1);
yna1= zeros(na,1);
unb = zeros(nb,1);
pnb = zeros(nbp,1);

%% Simulação

for T=1:N
    y_medido = -Denconv(1,2:na)*yna1(2:na,1)+Uconv*unb(1:nb,1)...
        +Pconv*pnb(1:nbp,1);
    yna = [zeros(1,1);yna(1:na-1,1)];
    yna(1) = y_medido;
    yna1 = [zeros(1,1);yna(1:na-1,1)];
    
    [f_u]=resp_livre(Ph, delta_u, y_medido, gd, salto);
    
    [fd]=resp_livre(Ph, delta_p, 0, dg, salto);
    
    wu=w(T+1:T+Ph,1);
    
    if T==1
        p_fut=p(T:T+Mh-1,1)-zeros(Mh,1);
    else
        p_fut=p(T:T+Mh-1,1)-p(T-1:T+Mh-2,1);
    end
    
    [u_fut]=alg_controle(K,wu,f_u,D,p_fut,fd);
    
    delta_u= [zeros(1,1);delta_u(1:salto,1)];
    
    delta_u(1)=u_fut(1);
    
    unb=[zeros(2,1);unb(2:nb-1,1)];
    u_atual=delta_u(1)+u_atual;
    unb(2)=u_atual;
    
    pnb=[zeros(2,1);pnb(2:nbp-1,1)];
    pnb(2)=p(T);
    
    delta_p=[zeros(1,1);delta_p(1:salto,1)];
    delta_p(1) = deltap1;
    deltap1 = p_fut(1);
    
    y_saida(T)=y_medido;
    U_passado(T)=u_atual;
end

%% Plot

figure(1)
subplot(2,1,1);
plot(tempo,wp)
hold on
plot(tempo,y_saida)
grid on
legend('T1ref(t)','T1(t)')
title('Temperatura do Tanque')
xlabel('t (s)')
ylabel('Temperatura')
subplot(2,1,2);
hold on
plot(tempo,pp)
title('Pertubação')
xlabel('t(s)')
ylabel('Abertura da válvula')
grid on

figure(2)
plot(tempo,U_passado)
title('Ação de Controle')
xlabel('t(s)')
ylabel('U')
grid on

