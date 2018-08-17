function solution =  MotorSpeedMechDyn(motor,voltage,frequency,prevRotorSpeed,momentInertia,loadTorque,varargin)

[~,y] = ode45(@(t,y) MechDEq(t,y,voltage,frequency,motor,momentInertia,loadTorque),[0 1],prevRotorSpeed*pi/30,varargin{:});

loadRotorSpeed = y(end)*30/pi;

EQSolutions = motor.EQSolutions(voltage,frequency,loadRotorSpeed);

solution = [loadRotorSpeed,EQSolutions(:,1),EQSolutions(:,5)];

end

function dwdt = MechDEq(~,rotorAngSpeed,voltage,frequency,motor,inerta,loadTorque)
    
    rotorSpeed = rotorAngSpeed*30/pi;

    motorTorque = motor.getTorque(voltage,frequency,rotorSpeed);
    
    dwdt = (motorTorque - loadTorque)/inerta;
end





