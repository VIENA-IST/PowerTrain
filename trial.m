clear all
close all
clc

% initialize variables
motor = InductionMachine;
momentInertia = 0.048087;
loadTorque = 20;

[t,y] = ode45(@(t,y) MechDEq(t,y,motor,momentInertia,loadTorque),[0 1],0);

plot(t,y*30/pi)

function dwdt = MechDEq(~,rotorAngSpeed,motor,inerta,loadTorque)
    
    rotorSpeed = rotorAngSpeed*30/pi;

    motorTorque = motor.getTorque(18,18,rotorSpeed);
    
    dwdt = (motorTorque - loadTorque)/inerta;
end



