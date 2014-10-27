function [] = md_sampler_tryptic(timepoints, data, options)

    parHistory = zeros(options.burnin + options.samples, numel(options.start));
    LLHistory = zeros(options.burnin + options.samples, 1);
    
    stepSize = 0.05*options.start;
   
    npars = numel(options.start) - 1;

    oldpars = options.start;
    newpars = oldpars;
    sim = simplified_model_tryptic(oldpars, timepoints);
    oldLL = LL_tryptic(sim, data, 0.05);

    accepted = 0;
    
    for i=1:options.burnin
        
        % propose new parameter values
        newpars = oldpars + stepSize.*randn(1,npars+1);
        
        % get simulated values
        sim = simplified_model_tryptic(newpars(1:npars), timepoints);
       
        % calculate LL
        newLL = LL_tryptic(sim, data, newpars(npars+1));
        
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
        newpars = oldpars + stepSize.*randn(1,npars+1);
        
        % get simulated values
        sim = simplified_model_tryptic(newpars(1:npars), timepoints);
       
        % calculate LL
        newLL = LL_tryptic(sim, data, newpars(npars+1));
        
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
    
    figure(103)
    gplotmatrix([parHistory, LLHistory])

    figure()
    plot(parHistory)
    
    figure
    plot(LLHistory)
    
    maxi = min(find(LLHistory == max(LLHistory)));
    parHistory(maxi,:)
    LLHistory(maxi)
    
    figure()
    plot(0:0.1:max(timepoints), simplified_model_tryptic(parHistory(maxi,1:npars), 0:0.1:max(timepoints)), 'r');
    hold on
    plot(timepoints, data, 'o')
    
    
    
end