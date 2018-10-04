function subplotpowercurrent(voltage,frequency,rotorSpeed,EQCSolutions)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

ax(1) = subplot(2,3,1:3);
ax(2) = subplot(2,3,4);
ax(3) = subplot(2,3,5);
ax(4) = subplot(2,3,6);

% Rotor speed vs. (POut,Eff)
yyaxis(ax(1),'left')
plot(ax(1),rotorSpeed,EQCSolutions(:,6)*1e-3)
ylabel(ax(1),'Output Power [kW]')

yyaxis(ax(1),'right')
plot(ax(1),rotorSpeed,EQCSolutions(:,end))
ylabel(ax(1),'eff [%]')

%Rotor speed vs. (Is,IM,Ir)
plot(ax(2),rotorSpeed,EQCSolutions(:,1))
ylabel(ax(2),'Stator Current [A]')
plot(ax(3),rotorSpeed,EQCSolutions(:,2))
ylabel(ax(3),'Mag. Current [A]')
plot(ax(4),rotorSpeed,EQCSolutions(:,3))
ylabel(ax(4),'Rotor Current [A]')

for aHandle = ax
    xlabel(aHandle,'Speed [RPM]')
    grid(aHandle,'minor')
end

titleString = sprintf('U = %g [V], F = %g [Hz], V/F = %0.2g',voltage/frequency, voltage/frequency);
title(ax(1),titleString)

end

