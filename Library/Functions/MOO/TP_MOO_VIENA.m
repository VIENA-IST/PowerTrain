function [y, cons] = TP_MOO_VIENA(x,motor,opt)

    %parametros:
    dTime = 1;
    time  = 0:dTime:30;
    timeLength = numel(time);

    %% variáveis de decisão;

    voltage   = x(1:timeLength); % composta
    frequency = x((timeLength+1):end);

    loadTorque = 40;

    dist = 0;
    distArray     = zeros(1,length(time));
    velocityArray = zeros(1,length(time));

    Energy=0;
    IStator = zeros(1,length(time));
    PStator = zeros(1,length(time));
    FLAG    = zeros(1,length(time));
    wRotor  = zeros(1,length(time));
    
    for i=1:length(time)

        [s, FLAG(i)]=motor.getRotorSpeed(voltage(i), frequency(i), wRotor(max(1,i-1))*30/pi, loadTorque,opt);
        if FLAG(i) == 0
            IStator(i)   = 0;
            PStator(i)   = 0;
            distArray(i) = dist;
        else
            wRotor(i)  = s(1);
            wWheels    = wRotor(i)*pi/30/8.6;
            IStator(i) = s(2);
            PStator(i) = s(3);
            
            
            velocityArray(i) = wWheels*0.25;
            
            dist = dist + velocityArray(i)*dTime;
            distArray(i) = dist;
        end
        Energy = Energy + PStator(i)*dTime;

    end
    %% Objectivos
    y = [0 0];
    cons = [0 0];

    % Distancia máxima
    y(1) = -distArray(end);

    % Minimizar energia consumida
    % y(2) = C;
    %y(2) = Energy;

    %% calculate the constraint violations
    Imax=120;
    
    c1 = Imax-max(IStator);
    c2 = nnz(~FLAG)*-1;
    
    if(c1<0)
        cons(1) = abs(c1);
    end
    
    if(c2<0)
        cons(2) = abs(c2);
    end
end



