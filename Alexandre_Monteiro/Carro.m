% CARRO
close all
clear all
clc

%% Alinea 1

% Tensao aplicada:
Us    = [76 121 121]; %[V]
freq  = [76 220 305]; %[Hz]
pp    = 2; % Dois pares de polos


% Dados do circuito equivalente:
Rs = 8.56e-3;
Ls = 0.067e-3;
Rr = 5.10e-3;
Lr = 0.06e-3;
Lm = 1.012e-3;
Rm = inf;


Zs = Rs+1j*(2*pi*freq)*Ls;
Zr = Rr+1j*(2*pi*freq)*Lr;
if Rm == inf
    Zm = 1j*(2*pi*freq)*Lm;
else
    Z1 = 1j*(2*pi*freq)*Lm;
    Z2 = Rm;
    Zm = (Z1*Z2)/(Z1+Z2);
end



for i=1:(length(Us)) % I1=S; I2=Mag.; I3 = R;
    X = strcat('f',num2str(freq(i)));
    U  = [Us(i); Us(i); 0]/(sqrt(3));
    ws = 2*pi*freq(i)/pp;
    wrr.(X) = 0:0.1:ws;
    wr = 0:0.1:ws;
    ss.(X) = (ws-wr)/ws;
    s = (ws-wr)/ws;
    
    for j=1:(length(s))
        R  = (Rr*(1-s(j)))/s(j);
        
        Z  = [Zs(i)  Zm(i)       0;
             Zs(i)    0   (Zr(i)+R);
             -1       1         1];

        I = Z\U; 
        
        
        Uout.(X)(j) = I(3)*R;
        iout.(X)(j) = I(3);
          im.(X)(j) = I(2);
          Um.(X)(j) = I(2)*Zm(i);
          is.(X)(j) = I(1);
        Pout.(X)(j) = 3*(abs(I(3)).^2)*R;
        
    end
            
end

%% Alinea 2
figure();
cor = NaN(length(Us),3);

for i=1:(length(Us))
    X = strcat('f',num2str(freq(i)));
    for j=1:(length(is.(X)))
        Pin.(X)(j,1) = 3*Us(i)*abs(is.(X)(j))*cos(angle(is.(X)(j)));
        ren.(X)(j,1) = Pout.(X)(j)/Pin.(X)(j);
        Tout.(X)(j,1) = Pout.(X)(j)/wrr.(X)(j);
    end
    cor(i,:) = rand(3,1);
    plot(wrr.(X),Tout.(X)','Color',cor(i,:));
    legend(X);
    hold on;
end
legend(num2str(freq'));
xlabel('Wr [rad/s]');
ylabel('T [N.m]');
hold off;


%% Alinea 3
clear all;
Ummax = 75; % [V] tensao maxima de magnetizacao para nao haver saturacao 

% Dados do circuito equivalente:
Rs = 9.679e-3;
Ls = 0.095e-3;
Rr = 10.754e-3;
Lr = 0.086e-3;
Lm = 1.461e-3;
Rm = inf;
pp    = 2;

freq = 5:5:350;
WR = 0:(2*pi*(max(freq))/pp);
TT = 0:450;
ren = NaN(length(WR),length(TT),length(freq));
RMSU = NaN(length(WR),length(TT),length(freq));
Umag = NaN(length(WR),length(TT),length(freq));
for z = 1:length(freq)
    ws = 2*pi*freq(z)/pp;    
    Zs = Rs+1j*(2*pi*freq(z))*Ls;
    Zr = Rr+1j*(2*pi*freq(z))*Lr;
    if Rm == inf
        Zm = 1j*(2*pi*freq(z))*Lm;
    else
        Z1 = 1j*(2*pi*freq(z))*Lm;
        Z2 = Rm;
        Zm = (Z1*Z2)/(Z1+Z2);
    end
    I = NaN(3,1);
    for i = 1:length(WR)
        for j = 1:length(TT)
            if WR(i) < ws
                P = (TT(j)*WR(i))/3;
                s = (ws-WR(i))/ws;
                R = Rr*((1-s)/s);
                I(3) = sqrt(P/R);
                I(2) = I(3)*(Zr+R)/Zm;
                Umag(i,j,z) = I(2)*Zm;
                I(1) = I(2)+I(3);
                u = I(1)*Zs+Umag(i,j,z);
                %if abs(Umag) < Ummax
                RMSU(i,j,z) = abs(u)*sqrt(3);
                ren(i,j,z) = (P/(abs(u)*abs(I(1))*cos(angle(u)-angle(I(1)))))*100;
                %end
            end
        end
    end
end
fbest = NaN(length(WR),length(TT));
ubest = NaN(length(WR),length(TT));
mbest = NaN(length(WR),length(TT));

for i = 1:length(WR)
    for j = 1:length(TT)
      [best,ind] = max(ren(i,j,:));
      fbest (i,j) = freq(ind);
      ubest (i,j) = RMSU(i,j,ind); 
      mbest (i,j) = Umag(i,j,ind); 
     end
end
     
[TT,WR] = meshgrid(TT,WR);

figure();
s1 = surf(TT,WR,ubest);
set(s1,'LineStyle','none');
ylabel('wr [rad/s]');
xlabel('T [N.m]');
zlabel('Us [V]');

figure();
s2 = surf(TT,WR,fbest);
set(s2,'LineStyle','none')
ylabel('wr [rad/s]');
xlabel('T [N.m]');
zlabel('Frequencia');

vf = (abs(ubest)./fbest);
ef = (abs(mbest)./fbest);

figure();
s3 = surf(TT,WR,vf);
set(s3,'LineStyle','none')
zlabel('V/f');
ylabel('wr [rad/s]');
xlabel('T [N.m]');

figure();
s4 = surf(TT,WR,ef);
set(s4,'LineStyle','none')
zlabel('E/f');
ylabel('wr [rad/s]');
xlabel('T [N.m]');


%% Alinea 4 - SIMULINK

% FIAT Seicento ELETTRA
Fiat = [ 1200;  % Weight [kg]
          400;  % Battery weight [kg]
           90;  % Autonomy [km]
          100;  % M speed [km/h]
            8;  % Acceleration (0-50km/h) [s]
        60.97;  % Wheel diameter [cm]
        150.8;  % Width [cm]
        144.5;  % Height [cm]
         2.18;  % Area [m^2]
         0.33 ]; % Drag coefficient
           
    


