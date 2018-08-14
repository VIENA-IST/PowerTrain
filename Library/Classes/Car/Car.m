classdef Car
    %CAR general mechanical parameters of a car.
    %   Parameters:
    %       Mass
    %       Front Area
    %       Tire Radius
    %       Wheel + Axis Inertia
    %       
    
    properties (Access = protected, Constant)
        Drag Coefficient = 1;
    end
    
    properties
        mass
        frontArea
        tireRadius
        totalInertia
        wheelInertia
        axisInertia
    end
    
    properties (Depedent)
        massEquivalent
    
    methods
        %Constructor
        function newCar = Car(mass,frontArea,tireRadius,wheelInertia,axisInertia)
            
        end
    end
    
end

