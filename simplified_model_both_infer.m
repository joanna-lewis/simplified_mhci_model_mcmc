timepoints = [0, 0.5, 1:8]; % times at which experimental observations were collected

options.npep = 4; % number of peptides followed

% Read in tryptic and epitope data.

tryptic_data = csvread('all_Kb_tryptic.csv',1,1);
Kb_data = csvread('all_Kb_epitope.csv',1,1);

% Set the number of burnin and sampling iterations for the sampler.
options.burnin = 1000;
options.samples = 1000;

% Initialise the sampler. These are good starting values for the example
% data:
options.start = [321.0953,  182.3905,  141.8622,  608.5598,... % A
                8.8215,    0.8412,    1.4392,    8.4441, ... % n
                1.1653,    0.1213,    0.4661,    0.7841, ... % ktr
                1588.9,    1344.5,    215.7,    2677.0, ... % u
                0.0019, ... % g1
                0.0109,    0.0220,    0.0466,    0.0143, ... % tryptic error
                38.3968,   30.9699,  299.7263,   30.2299 ... % epitope erroe
                ];


% Initialise the sampler step size. These are good starting values for the example
% data:
options.stepSize = [6.0920,    2.1632,    2.1497,    9.3343,... % A
                0.0658,    0.0086,    0.0319,    0.0371, ... % n
                0.0082,    0.0007,    0.0051,    0.0063, ... % ktr
                55.1929,   36.4525,    5.1300,   314.1086, ... % u
                0.0099, ... % g1
                0.0046,    0.0063,    0.0118,  0.0050, ... % tryptic error
                11.4384,   15.7938,  174.3392,   15.6905 ... % epitope error
                ];

% Now run the sampler:
gibbs_sampler(timepoints, tryptic_data, Kb_data, options)