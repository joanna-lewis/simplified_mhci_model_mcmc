timepoints = [0, 0.5, 1:8];

options.npep = 4; % no. peptides

tryptic_data = csvread('all_Kb_tryptic.csv',1,1);
%tryptic_data = tryptic_data(1:2,:)

Kb_data = csvread('all_Kb_epitope.csv',1,1);
%Kb_data = Kb_data(1:2,:)

options.burnin = 1000;
options.samples = 10000;

% options.start = 0.9*[3.7017e+02, 1.4098e+02, ... % A
%                 9.2331, 1.4801, ... % n
%                 1.0749, 5.0450e-01, ... % ktr
%                 1.0117e+02, 5.9327e+01, ... % u
%                 6.74, ... % g1
%                 0.01, 0.02,  ... % tryptic error
%                 5, 50, ... % epitope erroe
%                 ];
% 

% These are good starting for Kb
options.start = [317.0,    180.2,    164.2,    1024.4,... % A
                8.5188,    0.9097,    2.0417,   10.3491, ... % n
                1.1213,    0.1324,    0.5987,    0.9251, ... % ktr
                1211.8,    1074.6,    186.8,    2599.0, ... % u
                0.0453, ... % g1
                0.02, 0.02, 0.02,  0.2, ... % tryptic error
                50, 200, 500, 50 ... % epitope erroe
                ];


% These are good starting for Db
% options.start = 0.9*[3.7017e+02, 1.4098e+02, 1.6664e+02, ... % A
%                 9.2331, 1.4801, 1.8995, ... % n
%                 1.0749, 5.0450e-01, 4.1892e-01, ... % ktr
%                 1.0117e+02, 5.9327e+01, 7.3351e+01, ... % u
%                 6.74, ... % g1
%                 0.01, 0.02, 0.02,  ... % tryptic error
%                 5, 50, 10 ... % epitope erroe
%                 ];

% These are good step size for Db
options.stepSize = [9.2851,    2.9674,    1.9740,   10.4762,... % A
                0.0903,    0.0078,    0.0218,    0.0957, ... % n
                0.0154,    0.0012,    0.0052,    0.0087, ... % ktr
                42.3092,   45.4577,    7.0370,   82.7146, ... % u
                0.0161, ... % g1
                0.02, 0.02, 0.02,  0.2, ... % tryptic error
                50, 200, 500, 50 ... % epitope erroe
                ];

% These are good step size for Db
% options.stepSize = [3.6567, 3.6567, 3.6567, ...
%                 0.0384, 0.0384,  0.0384, ...
%                 0.0089, 0.0089, 0.0089, ...
%                 1.5388, 1.5388, 1.5388, ...
%                 0.1748, ...
%                 0.0026, 0.0052, 0.0052, ...
%                 1.2969, 12.9687, 10 ...
%                 ];


% options.stepSize = [0, 0.01*1.4098e+02, 0, ...
%                 0, 0, 0, ...
%                 0, 0, 0, ...
%                 0, 0, 0, ...
%                 0, ...
%                 0, 0, 0, ...
%                 0, 0, 0 ...
%                 ];
% 
            
% options.stepSize = 0.05*[0.05*3.7017e+02, 0.05*1.4098e+02, 0.05*1.6664e+02, ...
%                 0.05*9.2331, 0.05*1.4801, 0.05*1.8995, ...
%                 0.05*1.0749, 0.05*5.0450e-01, 0.05*4.1892e-01, ...
%                 0.05*1.0117e+02, 0.05*5.9327e+01, 0.05*7.3351e+01, ...
%                 0.05*6.74, ...
%                 0.01, 0.02, 0.02, ...
%                 5, 50, 10 ...
%                 ];
            
% options.parsToInfer

% for i=1:options.npep
%     
%     solstart = simplified_model(1, ...
%                     [options.start(i), ...
%                     options.start(options.npep + i), ...
%                     options.start(2*options.npep + i), ...
%                     options.start(3*options.npep + i), ...
%                     options.start(4*options.npep + 1) ...
%                     ], ...
%                 0:0.1:max(timepoints));
% 
%     figure(1)
%     subplot(1 ,options.npep, i)
%     plot(0:0.1:max(timepoints), solstart(1,:), 'r');
%     hold on
%     subplot(1 ,options.npep, i)
%     plot(timepoints, tryptic_data(i,:), 'o')
% 
%     figure(2)
%     subplot(1 ,options.npep, i)
%     plot(0:0.1:max(timepoints), solstart(2,:), 'r');
%     hold on
%     subplot(1 ,options.npep, i)
%     plot(timepoints, Kb_data(i,:), 'o')
% 
% end


gibbs_sampler(timepoints, tryptic_data, Kb_data, options)