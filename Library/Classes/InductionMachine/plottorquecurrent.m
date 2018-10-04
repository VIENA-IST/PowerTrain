function pHandle = plottorquecurrent(voltage,frequency,rotorSpeed,EQCSolutions)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[maxTorque,maxTorqueIndex] = max(EQCSolutions(:,4));
rotorSpeedMaximizer = rotorSpeed(maxTorqueIndex);

yyaxis left
pHandle(1:2) = plot(rotorSpeed',EQCSolutions(:,4),'b',rotorSpeedMaximizer,maxTorque,'*k');
ylabel('Torque [Nm]')

yyaxis right
pHandle(3) = plot(rotorSpeed',EQCSolutions(:,1),'r');
ylabel('Stator Current [A]')

grid minor
xlabel('speed [RPM]')
titleString = sprintf('U = %g [V], F = %g [Hz], V/F = %0.2g',voltage,frequency, voltage/frequency);
title(titleString)

end

