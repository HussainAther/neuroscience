addpath bdtoolkit-2018b/

% construct the model
sys = FHN();

% parameters
%sys.pardef(5).value = 1000;     % Idur = 1000
sys.pardef = bdSetValue(sys.pardef,'Idur',1000);

% time domain
sys.tspan = [0 200];

% solver options
sys.odeoption.RelTol = 1e-5;

% initial conditions
sys.vardef(1).value = 0;    % V=0
sys.vardef(2).value = 0;    % W=0
