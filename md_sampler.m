function [] = md_sampler(timepoints, tryptic_data, Kb_data, options)

    parHistory = zeros(options.burnin + options.samples, numel(options.start));
    LLHistory = zeros(options.burnin + options.samples, 1);
    
    stepSize = 0.05*options.start;
   
    npars = numel(options.start) - 2;

    oldpars = options.start;
    newpars = oldpars;
    sim = simplified_model(oldpars, timepoints);
    oldLL = LL(sim, tryptic_data, Kb_data, oldpars(npars+1:end));

    accepted = 0;
    
    for i=1:options.burnin
        
        % propose new parameter values
        newpars = oldpars + stepSize.*randn(1,npars+2);
        
        % get simulated values
        sim = simplified_model(newpars(1:npars), timepoints);
       
        % calculate LL
        newLL = LL(sim, tryptic_data, Kb_data, newpars(npars+1:end));
        
        % calculate acceptance ratio
        logRatio = newLL - oldLL;
        
        % accept or reject
        if logRatio > log(rand(1))
           
           oldpars = newpars; 
           oldLL = newLL;
           accepted = accepted + 1;

            
        end
        
        parHistory(i,:) = oldpars;
        LLHistory(i) = oldLL;
        
        if(mod(i,100) == 0)
                     
            accepted
            
            if(accepted > 50)
                stepSize = stepSize*1.1
            elseif(accepted < 20)
                stepSize = stepSize*0.9
            end
            
            accepted = 0;
            
        end
        
    end
    
    for i=1:options.samples
        
        % propose new parameter values
        newpars = oldpars + stepSize.*randn(1,npars+2);
        
        % get simulated values
        sim = simplified_model(newpars(1:npars), timepoints);
       
        % calculate LL
        newLL = LL(sim, tryptic_data, Kb_data, newpars(npars+1:end));
        
        % calculate acceptance ratio
        logRatio = newLL - oldLL;
        
        % accept or reject
        if logRatio > log(rand(1))
           
           oldpars = newpars; 
           oldLL = newLL;
            
        end
        
        parHistory(options.burnin+i,:) = oldpars;
        LLHistory(options.burnin+i) = oldLL;
        
    end
    
    figure(101)
    plot(parHistory)
    
    figure(102)
    plot(LLHistory)
    
    figure(103)
    gplotmatrix([parHistory, LLHistory])
    
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
    plot(0:0.1:max(timepoints), solmax(1,:), 'r');
    hold on
    plot(timepoints, tryptic_data, 'o')
    
    figure(2)
    plot(0:0.1:max(timepoints), solmax(2,:), 'r');
    hold on
    plot(timepoints, Kb_data, 'o')
    
end