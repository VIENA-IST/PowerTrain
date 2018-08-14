% Copyright (c) 2018 F. Ferreira da Silva
%
% Permission is hereby granted, free of charge, to any person obtaining a
% copy of this software and associated documentation files (the
% "Software"), to deal in the Software without restriction, including
% without limitation the rights to use, copy, modify, merge, publish,
% distribute, sublicense, and/or sell copies of the Software, and to permit
% persons to whom the Software is furnished to do so, subject to the
% following conditions:
%
% The above copyright notice and this permission notice shall be included
% in all copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
% NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
% DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
% OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
% USE OR OTHER DEALINGS IN THE SOFTWARE.

classdef InductionMachine
    %INDUCTIONMACHINE general description of Induction Machine.
    %   Gives nameplate information as well as derived information from
    %   nameplate data.
    %
    %   The electrical parameters are from the steinmetz equivalent
    %   circuit
    
    
%--------------------------------------------------------------------------
%-PROPERTIES---------------------------------------------------------------
%--------------------------------------------------------------------------

    properties (SetAccess = private)
        %NAMEPLATE DATA
        
        nameplateData
        
        %Structure containing:
        %   Frequency
        %   Voltage
        %   Current
        %   Power
        %   Speed
        %   PowerFactor
        %   Torque
        %   maxSpeed
        %   etc...
    end
        
    properties
        %STEINMETZ EQUIVALENT CIRCUIT PARAMETERS
        
        EQCircuit
        
        %Structure containing:
        %   stator Resistance
        %   stator LeakageResistance
        %   Magnetic Resistance
        %   Magnetic Inductance
        %   Rotor Resistance
        %   Rotor Leakage Inductance  
    end
    
    properties (Dependent)
        polePair
    end

%--------------------------------------------------------------------------
%-METHODS------------------------------------------------------------------
%--------------------------------------------------------------------------    
    
    methods

        %Constructor
        function newIM = InductionMachine(nameplateData,EQCircuit)
            switch nargin
                case 0
                    
                    newIM.nameplateData = ...
                        struct('Frequency',76,...
                               'Voltage',76,...
                               'Current',157,...
                               'Power',15e3,...
                               'Speed',2200,...
                               'PowerFactor',0.85,...
                               'Torque',65,...
                               'maxSpeed',10e3);

                    newIM.EQCircuit = ...
                              [8.56e-3,...    %Rs
                               0.06292e-3,...    %Ls
                               inf,...         %RM
                               1.01220e-3,...    %LM
                               5.10e-3,...   %Rr
                               0.06709e-3];      %Lr
                case 1
                    
                    %check if arguments are structures
                    flagStructND = isstruct(nameplateData);

                    assert(flagStructND,'Variable nameplateData not a structure.')

                    %check if nameplate data has at least the first 5 fields.
                    NDFieldNames = fieldnames(nameplateData);

                    flagFieldNamesND = ismember(NDFieldNames,...
                        {'Frequency','Voltage','Current','Power','Speed'});

                    assert(sum(flagFieldNamesND) == 5,...
                        'Structure nameplateData does not have one of the following fields: Frequency, Voltage, Current, Power or Speed.');


                    newIM.nameplateData = nameplateData;

                case 2
                    
                    %check if arguments are structures and array
                    flagStructND = isstruct(nameplateData);
                    flagArrayEC = isnumeric(EQCircuit);

                    assert(flagStructND,'Variable nameplateData not a structure.')
                    assert(flagArrayEC,'Variable EQCircuit not an array.')

                    %check if nameplate data has at least the first 5 fields.
                    NDFieldNames = fieldnames(nameplateData);

                    flagFieldNamesND = ismember(NDFieldNames,...
                        {'Frequency','Voltage','Current','Power','Speed'});

                    assert(sum(flagFieldNamesND) == 5,...
                        'Structure nameplateData does not have one of the following fields: Frequency, Voltage, Current, Power or Speed.');
                    
                    %check if EQCircuit is a 6-array
                    checkEQCDim = numel(EQCircuit);
                    assert(checkEQCDim == 6,...
                        'Missing %d element in variable EQCircuit',6-numel(EQCircuit))

                    newIM.nameplateData = nameplateData;
                    newIM.EQCircuit = EQCircuit;
                    
                otherwise
                    error('exceeded number of inputs')
                
            end
        end
        
%--------------------------------------------------------------------------
%-OVERLOAD-----------------------------------------------------------------
%--------------------------------------------------------------------------
        
        %DISPLAY
        function disp(newIM)
            
            fieldsIM = fieldnames(newIM);
            
            if any(ismember(fieldsIM,'nameplateData'))
                
                fprintf(1,'--------------------------\n')
                fprintf(1,'Nameplate Data\n')
                fprintf(1,'--------------------------\n')
                
                for field = fieldnames(newIM.nameplateData)'
                    fprintf(1,[field{1} ':% g\n'],newIM.nameplateData.(field{1}));
                end
                
            end
            
            if any(ismember(fieldsIM,'EQCircuit'))
                
                
                fprintf(1,'--------------------------\n')
                fprintf(1,'Equivalent Circuit Data\n')
                fprintf(1,'--------------------------\n')
                
                
                fieldName = {'Rs','Ls','RM','LM','Rr','Lr'};
                for field = 1:numel(fieldName)
                    fprintf(1,[fieldName{field} ':% g\n'],newIM.EQCircuit(field));
                end
                
                fprintf(1,'--------------------------\n\n')
                
            end
        end
                
%-------PLOT---------------------------------------------------------------
        function plot(newIM,varargin)
            
            %check how many input arguments are
            %set voltage, frequency and flag
            switch nargin
                case 1
                    voltage = newIM.nameplateData.Voltage;
                    frequency = newIM.nameplateData.Frequency;
                    FULL_FLAG = false;
                case 2
                    voltage = newIM.nameplateData.Voltage;
                    frequency = newIM.nameplateData.Frequency;
                    assert(strcmp(varargin{1},'full'),'flag not known. Did you mean: full')
                    FULL_FLAG = true;
                case 3
                    voltage = varargin{1};
                    frequency = varargin{2};
                    FULL_FLAG = false;
                case 4
                    voltage = varargin{1};
                    frequency = varargin{2};
                    assert(strcmp(varargin{3},'full'),'flag not known. Did you mean: full')
                    FULL_FLAG = true;
                otherwise
                    error('Too many input arguments')
            end
            
            %check if voltage and frequency are positive
            [voltage,frequency] = newIM.CheckVoltageFrequency(voltage,frequency);
            
            %calculate solutions for each rotor speed
            statorFreq = frequency/newIM.polePair;
            
            switch frequency
                case 0
                    rotorSpeed = 0;
                otherwise
                    rotorSpeed = (0:statorFreq/500:statorFreq)*60;
            end
            
            solutions = SweepSolutions(newIM,voltage,frequency,rotorSpeed);
            
            figure
            plottorquecurrent(voltage,frequency,rotorSpeed,solutions);
            
            if FULL_FLAG == true
                figure;
                subplotpowercurrent(voltage,frequency,rotorSpeed,solutions);
            end
            
        end
%-------END-PLOT-----------------------------------------------------------

%--------------------------------------------------------------------------
%-SET METHODS--------------------------------------------------------------
%--------------------------------------------------------------------------

        
        function newIM = set.EQCircuit(newIM,EQCircuit)
            
            %check if EQCircuit is a strucutre
            flagStructEC = isnumeric(EQCircuit);
            assert(flagStructEC,'Variable EQCircuit not an array.')
            
            %check if EQCircuit is a 6-array
            checkEQCDim = numel(EQCircuit);
            assert(checkEQCDim == 6,'Missing %d element in variable EQCircuit',6-numel(EQCircuit))
            
            newIM.EQCircuit = EQCircuit;
        end
        
%--------------------------------------------------------------------------
%-GET METHODS--------------------------------------------------------------
%--------------------------------------------------------------------------
        
        %POLEPAIR
        function polePair = get.polePair(newIM)
            nomFrequencyRPM = newIM.nameplateData.Frequency*60;
            nomSpeed = newIM.nameplateData.Speed;
            polePair = round(nomFrequencyRPM/nomSpeed);
        end
        
%--------------------------------------------------------------------------
%-PUBLIC METHODS-----------------------------------------------------------
%--------------------------------------------------------------------------


%-------EQSOLUTIONS--------------------------------------------------------

        %Calculations using equivalent circuit (ADAPTED FROM ALEXANDRE MONTEIRO'S WORK)
        function solution = EQSolutions(newIM,voltage,frequency,rotorSpeed)

            %check if EQCircuit exists
            fieldsIM = fieldnames(newIM);
            assert(any(ismember(fieldsIM,'EQCircuit')),...
                'equivalent circuit parameters are needed (missing EQCircuit structure')
            
            %check if operationValues are not empty and are numeric
            assert(all(isnumeric([voltage frequency rotorSpeed])),...
                'Operation Values are not numeric.')
            
            %if frequency is zero, the solutions are 0
            if frequency == 0
                solution = zeros(1,7);
                return
            end
            
            %make calculations
            angFreq = 2*pi*frequency;
            statorFreq = frequency*60/newIM.polePair;
            slip = (statorFreq - rotorSpeed)/(statorFreq);
            
            eqvC = newIM.EQCircuit;
            % eqvC(1) == Rs 
            % eqvC(2) == Ls
            % eqvC(3) == RM
            % eqvC(4) == LM
            % eqvC(5) == Rr
            % eqvC(6) == Lr
            
            %Compute Z Matrix of the EQCircuit
            % U = ZI
            % U = [Us Us  0]
            % I = [Is IM Ir]
            % Z = [Zs ZM  0
            %      Zs  0 Zr
            %      -1  1  1]
            
            Zs = eqvC(1)+1j*angFreq*eqvC(2);
            
            %check if RM exists or not
            if eqvC(3) == inf
                Zm = 1j*angFreq*eqvC(4);
            else
                Z1 = 1j*angFreq*eqvC(4);
                Z2 = eqvC(3);
                Zm = (Z1*Z2)/(Z1+Z2);
            end
            
            switch slip
                % no load condition
                case 0
                    
                    U = voltage/sqrt(3);
                    Z = Zs + Zm;
                    I = [Z\U;Z\U;0];
                    
                    outputPower = 0;
                    Torque= 0;
                    
                otherwise
                    
                    U = [voltage voltage 0]/sqrt(3);
                    Zr = eqvC(5)/slip+1j*angFreq*eqvC(6);
                    Z = [Zs     Zm      0;
                         Zs     0       Zr;
                         -1     1       1];
                    I = Z\U';
                    
                    outputPower = 3*(abs(I(3))^2)*eqvC(5)*(1-slip)/slip;
                    Torque = 3*(abs(I(3))^2)*eqvC(5)/(slip*(angFreq/newIM.polePair));
            end
            
            inputPower = real(3*U(1)*conj(I(1)));
            efficiency = outputPower/inputPower;
            
            solution = [abs(I(1)),...    %Is
                        abs(I(2)),...    %IM
                        abs(I(3)),...    %Ir
                        Torque,...       %T
                        inputPower,...   %Pin
                        outputPower,...  %Pout
                        efficiency];     %eff            
        end
%-------END-EQSOLUTIONS----------------------------------------------------
        

%-------GETROTORSPEED------------------------------------------------------
        function [solution, STABLE_FLAG] = getRotorSpeed(newIM,voltage,frequency,prevRotorSpeed,loadTorque,varargin) %need to check negative values of rotor speed!!
            
            %check if voltage and frequency are positive
            [voltage,frequency] = newIM.CheckVoltageFrequency(voltage,frequency);
            
            %Pre-conditioned solution
            solution = [];
            STABLE_FLAG = 0;
            
            %check if there is a non-zero frequency
            if frequency == 0
                return
            end
            
            %get starting torque and rotorspeed maximizer
            [~, maximRotorSpeed] = newIM.getMaxTorque(voltage,frequency);
            
            %convert to slip and make slip condition
            maximSlip = newIM.RotorSpeed2slip(frequency,maximRotorSpeed);
            prevSlip = newIM.RotorSpeed2slip(frequency,prevRotorSpeed);
            
            %initialize function handles to fsolve
            torqueHandle = newIM.TorqueFunctionHandle(voltage,frequency);
            objFunction = @(s) abs(torqueHandle(s) - loadTorque);
            
            
            
            
            %if in stable region, check intersection
            slipCondition = prevSlip < maximSlip;
            torqueDifference = torqueHandle(prevSlip) - loadTorque;
            if slipCondition || (torqueDifference > 0)
                loadSlip = fminbnd(objFunction,0,maximSlip,varargin{:});
                loadRotorSpeed = newIM.Slip2RotorSpeed(frequency,loadSlip);
            else
                return
            end
            
            if ~isnan(loadRotorSpeed)
                EQSolutions = newIM.EQSolutions(voltage,frequency,loadRotorSpeed);
                solution = [loadRotorSpeed, EQSolutions(:,1), EQSolutions(:,5)];
                STABLE_FLAG = 1;
            end
        end
%-------END-GETROTORSPEED--------------------------------------------------

%-------GETSTALLTORQUE-----------------------------------------------------
        function startingTorque = getStartingTorque(newIM,voltage,frequency)
            
            %check if voltage and frequency are positive
            [voltage,frequency] = newIM.CheckVoltageFrequency(voltage,frequency);
            
            %initialize variables
            t = 2*pi;
            we = t*frequency;
            pp = newIM.polePair;
            ws = t*frequency/pp;
            
            eqc = newIM.EQCircuit;
            Rr = eqc(5);
            Xr = we*eqc(6);
            
            [VTH , ZTH] = newIM.getThevenin(voltage,frequency);
            
            Ir = abs(VTH/(ZTH + 1j*Xr + Rr));
            
            startingTorque = 3/ws * Ir^2 * Rr;
        end
%-------END-GETSTALLTORQUE-------------------------------------------------

%-------GETMAXTORQUE-------------------------------------------------------
        function [maxTorque, maximRotorSpeed] = getMaxTorque(newIM,voltage,frequency)
            
            %check if voltage and frequency are positive
            [voltage,frequency] = newIM.CheckVoltageFrequency(voltage,frequency);
            
            %initialize variables
            t = 2*pi;
            we = t*frequency;
            pp = newIM.polePair;
            ws = t*frequency/pp;
            
            eqc = newIM.EQCircuit;
            Rr = eqc(5);
            Xr = we*eqc(6);
            
            [VTH , ZTH] = newIM.getThevenin(voltage,frequency);
            
            RTH = real(ZTH);
            XTH = imag(ZTH);
            
            maxTorque = 3/(2*ws) * abs(VTH)^2 * 1/(RTH + sqrt(RTH^2 + (XTH + Xr)^2));
            
            maximSlip = Rr * 1/(RTH + sqrt(RTH^2 + (XTH + Xr)^2));
            
            maximRotorSpeed = newIM.Slip2RotorSpeed(frequency,maximSlip);
        end
%-------END-GETMAXTORQUE---------------------------------------------------

    end

%--------------------------------------------------------------------------
%-PRIVATE METHODS----------------------------------------------------------
%--------------------------------------------------------------------------
    
    methods (Access = private)
        
%-------CHECKIFPOSITIVE----------------------------------------------------
        function bool = CheckIfPositive(~,variable)
            if variable >=0
                bool = true;
            else
                bool = false;
            end
        end
%-------END-CHECKIFPOSITIVE------------------------------------------------

%-------CHECKVOLTAGEFREQUENCY----------------------------------------------
        function [v,f] = CheckVoltageFrequency(newIM,voltage,frequency)
            
            v = voltage;
            f = frequency;
            
            if ~newIM.CheckIfPositive(frequency)
                warning('frequency negative. changing it to positive')
                f = -frequency;
            end
            if ~newIM.CheckIfPositive(voltage)
                warning('voltage negative. changing it to positive')
                v = -voltage;
            end
        end
%-------END-CHECKIFPOSITIVE------------------------------------------------

%-------SWEEPSOLUTIONS-----------------------------------------------------
        function solutions = SweepSolutions(newIM,voltage,frequency,rotorSpeed)
                solutions = NaN(numel(rotorSpeed),7);
            for k = 1:numel(rotorSpeed)
                solutions(k,:) = EQSolutions(newIM,voltage,frequency,rotorSpeed(k));
            end
        end
%-------END-SWEEPSOLUTIONS-------------------------------------------------

%-------SLIP2ROTORSPEED----------------------------------------------------
        function rotorSpeed = Slip2RotorSpeed(newIM,frequency,slip)
                pp = newIM.polePair;
                
                statorSpeed = frequency*60/pp;
                
                rotorSpeed = (1-slip)*statorSpeed;
        end
%-------END-SLIP2ROTORSPEED------------------------------------------------

%-------SLIP2ROTORSPEED----------------------------------------------------
        function slip = RotorSpeed2slip(newIM,frequency,RotorSpeed)
                pp = newIM.polePair;
                
                statorSpeed = frequency*60/pp;
                
                slip = 1 - RotorSpeed/statorSpeed;
        end
%-------END-SLIP2ROTORSPEED------------------------------------------------
    
%-------GETTHEVENIN--------------------------------------------------------
        function [VTH,ZTH] = getThevenin(newIM,voltage,frequency)
            
            v = voltage/sqrt(3);
            we = 2*pi*frequency;
            
            eqc = newIM.EQCircuit;           
            Rs = eqc(1);
            Xs = we*eqc(2);
            XM = we*eqc(4);
            
            VTH = v*1j*XM/(Rs + 1j*(Xs + XM));
            ZTH = 1j*XM*(Rs + 1j*Xs)/(Rs + 1j*(Xs + XM));
        end
%-------END-GETTHEVENIN----------------------------------------------------

    
%-------GETTHEVENIN--------------------------------------------------------
        function fHandle = TorqueFunctionHandle(newIM,voltage,frequency)
        
            %check if voltage and frequency are positive
            [voltage,frequency] = newIM.CheckVoltageFrequency(voltage,frequency);
            
            %initialize variables
            t = 2*pi;
            we = t*frequency;
            pp = newIM.polePair;
            ws = t*frequency/pp;
            
            eqc = newIM.EQCircuit;
            Rr = eqc(5);
            Xr = we*eqc(6);
            
            [VTH , ZTH] = newIM.getThevenin(voltage,frequency);
            
            Ir = @(s) abs(VTH/(ZTH + 1j*Xr + Rr/s));
            
            fHandle = @(s) 3/ws * (Ir(s))^2 * Rr/s;
        end
%-------END-GETTHEVENIN----------------------------------------------------

    end    
end

