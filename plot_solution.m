clear all
close all
clc

A=textread('Population_last.txt');
A=sortrows(A,63);


V_v=A(:,1:31);
f_v=A(:,32:62);
x_max_v=A(:,63);
E_v=A(:,64);
cons1_v=A(:,65);
cons2_v=A(:,66);

subplot(4,1,1);
plot(V_v);
subplot(4,1,2);
plot(f_v);
subplot(4,1,3);
plot(x_max_v);
subplot(4,1,4);
plot(E_v);

figure(2)
subplot(2,1,1);
plot(cons1_v);
subplot(2,1,2);
plot(cons2_v);


%%
%parametros:
dt=1;
t=0:dt:30;

%% variáveis de decisão;

% B=[122.779	163.955	92.2746	112.668	116.868	108.607	70.7539	58.2586	64.2967	40.7842	158.7	96.87	135.994	128.224	155.598	123.516	109.233	131.261	103.067	138.76	82.7746	53.6451	28.5687	15.7411	62.3245	137.828	86.3155	126.599	127.655	156.739	92.6681	49.8826	7.88974	11.7329	34.4649	5.59718	26.7448	21.2399	42.3924	3.53549	20.4144	36.694	37.3206	30.8977	16.8406	22.8612	29.7185	20.4757	7.79695	38.7375	8.6266	22.6777	3.46187	12.3628	12.1599	29.6639	34.4874	28.4718	24.2759	10.4014	23.2664	50	-234.27	345045	0	];
% V1=B(1:31);
% f1=B(32:62);

sol=1;
V=V_v(sol,:)'; % composta
f=f_v(sol,:)';

figure(3)
subplot(1,2,1)
plot(V)
subplot(1,2,2)
plot(f)

% import motor library
motor=InductionMachine;
%

% V=ones(1,length(t))*19*sqrt(3);
% f=ones(1,length(t))*33;

T_load=40;

x=0;
x_v=zeros(1,length(t));
v_v=zeros(1,length(t));

E=0;
Is=zeros(1,length(t));
Ps=zeros(1,length(t));
flag=zeros(1,length(t));
wr=zeros(1,length(t));
v_f=zeros(1,length(t));
for i=1:length(t)
    v_f(i)=V(i)/f(i);
    [s, flag(i)]=motor.getRotorSpeed(V(i), f(i), wr(max(1,i-1))*60/(2*pi), T_load);
    if flag(i)==0
        x=x;
        Is(i)=0;
        Ps(i)=0;
        x_v(i)=x;
    else
        wr(i)=s(1);
        wr_wheels=(s(1)*2*pi/60)/8.6;
        Is(i)=s(2);
        Ps(i)=s(3);
        v=wr_wheels*0.25;
        v_v(i)=v;
        x=x+v*dt;
        x_v(i)=x;
    end
    E=E+Ps(i)*dt;
    
end

figure(3)
subplot(1,4,1)
plot(t,x_v);
subplot(1,4,2)
plot(t,V);
subplot(1,4,3)
plot(t,Is);
subplot(1,4,4)
plot(t,flag);



