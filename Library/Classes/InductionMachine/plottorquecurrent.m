function pHandle = plottorquecurrent(voltage,frequency,rotorSpeed,EQCSolutions)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

maxTorqueIndex = find(diff(EQCSolutions(:,4)) >= 0);

if isempty(maxTorqueIndex)
    rotorSpeedMaximizer = 0;
    maxTorque = EQCSolutions(1,4);
else
    rotorSpeedMaximizer = rotorSpeed(maxTorqueIndex(end)+1);
    maxTorque = EQCSolutions(maxTorqueIndex(end)+1,4);
end

yyaxis left
pHandle(1:2) = plot(rotorSpeed',EQCSolutions(:,4),'b',rotorSpeedMaximizer,maxTorque,'*k');
ylabel('Torque [Nm]')

yyaxis right
pHandle(3) = plot(rotorSpeed',EQCSolutions(:,1),'r');
ylabel('Stator Current [A]')

grid minor

xlabel('speed [RPM]')

titleString = sprintf('Us = %g [V], F = %g [Hz], V/F = %0.2g',voltage,frequency, voltage/frequency);
title(titleString)


end

