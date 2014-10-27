function [] = gibbs_sampler(timepoints, tryptic_data, Kb_data, options)

    parHistory = zeros(options.burnin + options.samples, numel(options.start));
    LLHistory = zeros(options.burnin + options.samples, 1);
    
    stepSize = 0.05*options.start;
   
    npars = numel(options.start) - 2;

    oldpars = options.start;
    newpars = oldpars;
    sim = simplified_model(oldpars, timepoints);
    oldLL = LL(sim, tryptic_data, Kb_data, oldpars(npars+1:end));

    accepted = zeros(size(options.start));
    
    for i=1:options.burnin
        
        %i
        
        for j=1:numel(options.start)
            
            newpars = oldpars;
        
            % propose new parameter values
            newpars(j) = oldpars(j) + stepSize(j)*randn(1);
            % try again if you picked a parameter value < 0
            while min(newpars) < 0
                %newpars
                newpars(j) = oldpars(j) + stepSize(j)*randn(1);
            end

            % get simulated values
            sim = simplified_model(newpars(1:npars), timepoints);

            % calculate LL
            newLL = LL(sim, tryptic_data, Kb_data, newpars(npars+1:end));

            % calculate acceptance ratio
            logRatio = newLL - oldLL;

            % accept or reject
            if logRatio > log(rand(1))

               oldpars(j) = newpars(j); 
               oldLL = newLL;
               accepted(j) = accepted(j) + 1;
               
            end
        
        end
        
        parHistory(i,:) = oldpars;
        LLHistory(i) = oldLL;
        
        if(mod(i,100) == 0)
                     
            accepted
            
            stepSize(accepted > 50) = stepSize(accepted > 50)*1.1;
            
            stepSize(accepted < 20) = stepSize(accepted < 20)*0.9;
            
            accepted = zeros(size(options.start));
            
        end
        
        
    end
    
    for i=1:options.samples
        
        for j=1:numel(options.start)
            
            newpars = oldpars;

            % propose new parameter values
            newpars(j) = oldpars(j) + stepSize(j)*randn(1);
            % try again if you picked a parameter value < 0
            while min(newpars) < 0
                %newpars
                newpars(j) = oldpars(j) + stepSize(j)*randn(1);
            end
            
            % try again if you picked a parameter < 0
            while min(newpars < 0)
                newpars(j) = oldpars(j) + stepSize(j)*randn(1);
            end

            % get simulated values
            sim = simplified_model(newpars(1:npars), timepoints);

            % calculate LL
            newLL = LL(sim, tryptic_data, Kb_data, newpars(npars+1:end));

            % calculate acceptance ratio
            logRatio = newLL - oldLL;

            % accept or reject
            if logRatio > log(rand(1))

               oldpars(j) = newpars(j); 
               oldLL = newLL;

            end
            
        end
        
        parHistory(options.burnin+i,:) = oldpars;
        LLHistory(options.burnin+i) = oldLL;
        
    end
    
    figure(101)
    plot(parHistory(options.burnin+1:end,:))
    
    figure(102)
    plot(LLHistory(options.burnin+1:end))
    
    figure(103)
    gplotmatrix([parHistory(options.burnin+1:end,:), LLHistory(options.burnin+1:end)])
    
    maxi = min(find(LLHistory == max(LLHistory)));
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
                
    solmax = simplified_model(maxpars, 0:0.1:max(timepoints));
    
    figure(1)
    plot(0:0.1:max(timepoints), solmax(1,:), 'r --');
    hold on
    plot(timepoints, tryptic_data, 'o')
    
    figure(2)
    plot(0:0.1:max(timepoints), solmax(2,:), 'r --');
    hold on
    plot(timepoints, Kb_data, 'o')
    
end