
figure(101)
plot(parHistory(options.burnin+1:end,:))

figure(102)
plot(LLHistory(options.burnin+1:end))

figure(103)
gplotmatrix([parHistory(options.burnin+1:end,:), LLHistory(options.burnin+1:end)])


timepoints = [0, 0.5, 1:8];


maxi = min(find(LLHistory == max(LLHistory)));

maxpars = parHistory(maxi,:);

solmax = simplified_model(options.npep, maxpars, 0:0.1:max(timepoints));

h1 = figure;
set(h1, 'Units', 'inches', 'Position', [0,4.5, 4.5* options.npep, 4.5], 'PaperSize', [4.5* options.npep, 4.5])

h2 = figure;
set(h2, 'Units', 'inches', 'Position', [0,0, 4.5* options.npep, 4.5], 'PaperSize', [4.5* options.npep, 4.5])

% add some noise to simulations from sampled parameters, to give a
% simulated dataset
sims = simHistory(:,:,options.burnin+1:end);
sims = sims + repmat((maxpars((4*options.npep+2) : end))', [1, numel(timepoints), options.samples]).*randn(size(sims));


for i=1:options.npep

figure(h1)
subplot(1 ,options.npep, i)
fill([0, timepoints, fliplr(timepoints), 0], [0, quantile(sims(i,:,:), 0.025, 3), fliplr(quantile(sims(i,:,:), 0.975, 3)), 0], ...
    [1, 0.78, 0.80], ...
    'EdgeColor', 'none' ...
    )
hold on
plot(0:0.1:max(timepoints), solmax(i,:), 'r');
plot(timepoints, tryptic_data(i,:), 'o')

ylim([0, max(max(max(sims(1:3, :, :))))])

if (i==1)
ylabel('Epitope level (% maximum)')
end

if (i==2)
xlabel('Time post infection (hours)')
end

[0, timepoints, fliplr(timepoints), 0]
[0, quantile(sims(options.npep+i,:,:), 0.025, 3), fliplr(quantile(sims(options.npep+i,:,:), 0.975, 3)), 0]

figure(h2)
subplot(1 ,options.npep, i)
fill([0, timepoints, fliplr(timepoints), 0], [0, quantile(sims(options.npep+i,:,:), 0.025, 3), fliplr(quantile(sims(options.npep+i,:,:), 0.975, 3)), 0], ...
    [1, 0.78, 0.80], ...
    'EdgeColor', 'none' ...
    )
hold on
plot(0:0.1:max(timepoints), solmax(options.npep+i,:), 'r');
plot(timepoints, Kb_data(i,:), 'o')

if (i==1)
ylabel('Tryptic level (copies per cell)')
end

if (i==2)
xlabel('Time post infection (hours)')
end
    
end