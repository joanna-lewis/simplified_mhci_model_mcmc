function [] = gibbs_sampler(timepoints, tryptic_data, Kb_data, options)

    parHistory = zeros(options.burnin + options.samples, numel(options.start));
    LLHistory = zeros(options.burnin + options.samples, 1);
    simHistory = zeros(2*size(tryptic_data,1), numel(timepoints), options.burnin + options.samples);
    
    stepSize = options.stepSize;
   
    npars = numel(options.start) - 2*options.npep;

    oldpars = options.start;
    newpars = oldpars;
    
    oldsim = simplified_model(options.npep, oldpars(1:npars), timepoints)
    newsim = oldsim;
    
    oldLL = LL(oldsim, tryptic_data, Kb_data, oldpars(npars+1:end));

    accepted = zeros(size(options.start));
    
    for i=1:options.burnin
        
        %i
        
        for j=1:numel(options.start)
            
            %j
            
            newpars = oldpars;
        
            % propose new parameter values
            newpars(j) = oldpars(j) + stepSize(j)*randn(1);
            
%             oldpars(j)
%             newpars(j)
            
            % try again if you picked a parameter value < 0
            while min(newpars) < 0
                %newpars
                newpars(j) = oldpars(j) + stepSize(j)*randn(1);
            end

            if (j <= 4*options.npep) % if the parameter a peptide-specific model parameter

                pnum = mod(j, options.npep);
                if pnum == 0
                    pnum = options.npep;
                end
                
                newsim([pnum, options.npep + pnum],:) = simplified_model(1, newpars([pnum, options.npep+pnum, 2*options.npep+pnum, 3*options.npep+pnum, 4*options.npep+1]), timepoints);

            elseif (j == 4*options.npep + 1) % if it's g1
                
                % simulate all peptides
                newsim = simplified_model(options.npep, newpars(1:npars), timepoints);
                
            end
                
            % calculate LL
            %oldLL
            newLL = LL(newsim, tryptic_data, Kb_data, newpars(npars+1:end));

            % calculate acceptance ratio
            logRatio = newLL - oldLL;

            % accept or reject
            if logRatio > log(rand(1))

               oldpars(j) = newpars(j); 
               oldLL = newLL;
               oldsim = newsim;
               accepted(j) = accepted(j) + 1;
               
            end
        
        end
        
        parHistory(i,:) = oldpars;
        LLHistory(i) = oldLL;
        simHistory(:,:,i) = oldsim;
        
        if(mod(i,100) == 0)
            
            i
            accepted
            
            stepSize(accepted > 50) = stepSize(accepted > 50)*1.1;
            
            stepSize(accepted < 20) = stepSize(accepted < 20)*0.9;
            
            accepted = zeros(size(options.start));
            
        end
        
        
    end
    
    
    for i=1:options.samples
        
        %i
        
        for j=1:numel(options.start)
            
            %j
            
            newpars = oldpars;
        
            % propose new parameter values
            newpars(j) = oldpars(j) + stepSize(j)*randn(1);
            % try again if you picked a parameter value < 0
            while min(newpars) < 0
                %newpars
                newpars(j) = oldpars(j) + stepSize(j)*randn(1);
            end

            if (j <= 4*options.npep) % if the parameter a peptide-specific model parameter

                pnum = mod(j, options.npep);
                if pnum == 0
                    pnum = options.npep;
                end
                
            newsim([pnum, options.npep + pnum],:) = simplified_model(1, newpars([pnum, options.npep+pnum, 2*options.npep+pnum, 3*options.npep+pnum, 4*options.npep+1]), timepoints);

            elseif (j == 4*options.npep + 1) % if it's g1
                
                % simulate all peptides
                newsim = simplified_model(options.npep, newpars(1:npars), timepoints);
                
            end
                
            % calculate LL
            newLL = LL(newsim, tryptic_data, Kb_data, newpars(npars+1:end));

            % calculate acceptance ratio
            logRatio = newLL - oldLL;

            % accept or reject
            if logRatio > log(rand(1))

               oldpars(j) = newpars(j); 
               oldLL = newLL;
               oldsim = newsim;
               accepted(j) = accepted(j) + 1;
               
            end
        
        end
        
        parHistory(options.burnin + i,:) = oldpars;
        LLHistory(options.burnin + i) = oldLL;
        simHistory(:,:,options.burnin + i) = oldsim;
         
        if(mod(i,100) == 0)
            
            i
%             accepted
%             
%             stepSize(accepted > 50) = stepSize(accepted > 50)*1.1;
%             
%             stepSize(accepted < 20) = stepSize(accepted < 20)*0.9;
%             
%             accepted = zeros(size(options.start));
            
        end
        
    end
    
    % Save data?
    
    if strcmp(input('Save results? (y/n)','s') , 'y');
                
        save(['./Results/' num2str(fix(clock))], ...
            'i', ...
            'tryptic_data', ...
            'Kb_data', ...
            'options', ...
            'stepSize', ...
            'parHistory', ...
            'LLHistory', ...
            'simHistory' ...
            );
        
    end
    
    figure(101)
    plot(parHistory(options.burnin+1:end,:))
    
    figure(102)
    plot(LLHistory(options.burnin+1:end))
    
    figure(103)
    gplotmatrix([parHistory(options.burnin+1:end,:), LLHistory(options.burnin+1:end)])
    
    maxi = min(find(LLHistory == max(LLHistory)));
    
    %stepSize
    parHistory(maxi,:)
    LLHistory(maxi)
    
%     maxpars = [2.7, 3.1, 0.6, ... %A, n, ktr
%                 3600*0.1, ...            % e
%                 0.02, 20, ...   % g0, g1
%                 3600*10^-4, ...          % u
%                 0.1 ...           % D
%                 ];
    maxpars = parHistory(maxi,:);
%     maxpars(6) = parHistory(maxi,2);
                
    solmax = simplified_model(options.npep, maxpars, 0:0.1:max(timepoints));
    
    % add some noise to simulations from sampled parameters, to give a
    % simulated dataset
    sims = simHistory(:,:,options.burnin+1:end);
    sims = sims + repmat((maxpars((4*options.npep+2) : end))', [1, numel(timepoints), options.samples]).*randn(size(sims));

    for i=1:options.npep
    
    figure(1)
    subplot(1 ,options.npep, i)
    fill([timepoints, fliplr(timepoints)], [quantile(sims(i,:,:), 0.025, 3), fliplr(quantile(sims(i,:,:), 0.975, 3))], ...
        [1, 0.78, 0.80], ...
        'EdgeColor', 'none' ...
        )
    hold on
    plot(0:0.1:max(timepoints), solmax(i,:), 'r');
    plot(timepoints, tryptic_data(i,:), 'o')
    
    figure(2)
    subplot(1 ,options.npep, i)
    fill([timepoints, fliplr(timepoints)], [quantile(sims(options.npep+i,:,:), 0.025, 3), fliplr(quantile(sims(options.npep+i,:,:), 0.975, 3))], ...
        [1, 0.78, 0.80], ...
        'EdgeColor', 'none' ...
        )
    hold on
    plot(0:0.1:max(timepoints), solmax(options.npep+i,:), 'r');
    plot(timepoints, Kb_data(i,:), 'o')
    
    end