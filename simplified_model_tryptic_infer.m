timepoints = [0.5, 1:8];

data = csvread('B8_tryptic.csv',1,1);

options.burnin = 1000;
options.samples = 10000;
options.start = [2.7, 3.1, 0.6, 0.05];

md_sampler_tryptic(timepoints, data, options)