%*************************************************************************
% Test Problem: 'MOO_VF'
% Description: Optimization of scalar command V/f
%   (1)constrained
%
% Reference : [1] Deb K, Pratap A, Agarwal S, et al. A fast and elitist 
%   multiobjective genetic algorithm NSGA-II[J]. Evolutionary Computation. 
%   2002, 6(2): 182-197.
%*************************************************************************

%% Initialize variables and fminbnd options


motor = InductionMachine;
loadTorque = 40;
Inertia = 1.48;

opt = odeset('RelTol',1e-3,'AbsTol',1e-3);

%% NSGA Config


options = nsgaopt();                    % create default options structure
options.popsize = 500;                   % populaion size
options.maxGen  = 100;                  % max generation

options.numObj = 1;                     % number of objectives
options.numVar = 2*31;                  % number of design variables
options.numCons = 1;                    % number of constraints

options.vartype = ones(1,options.numVar);

% voltage frequency
options.lb = [ones(1,31)*0.1 ones(1,31)*0.1];                  % lower bound of x
options.ub = [ones(1,31)*200 ones(1,31)*300];                  % upper bound of x

% parallel computing
options.useParallel = 'yes';
options.poolsize = 3;

options.objfun = @TP_MOO_VIENA;         % objective function handle
options.plotInterval = 5;               % interval between two calls of "plotnsga". 

%% NSGA Run


result = nsga2(options,motor,Inertia,loadTorque,opt);          % begin the optimization!


