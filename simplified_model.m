function[sol] = simplified_model(npep, parameters, timepoints)

    % simplified model of peptide presentation
    
    % npep is the number of peptides
    
    % ode solver options
    options = odeset('RelTol', 10^-10, 'AbsTol', 10^-10);

    
    % set up vector for solution
    sol = zeros(2*npep, numel(timepoints));
    
    for i=1:npep
    
        % allocate parameters
        A = parameters(i);
        n = parameters(npep + i);
        ktr = parameters(2*npep + i);
        
        u = parameters(3*npep + i); %5.9327e+01; %3600*10^-4;

        e = 1e4;

        g0 = 1;
        g1 = parameters(4*npep + 1); %6.74;

        % peptide levels
        sol(i,:) = 0.01*peptide(timepoints);
        
        [~, sol(npep+i,:)] = ode15s(@ydot, timepoints, 0, options);
        %sol(2,:) = sol(2,:);
        
    end
    
    % function for calculating derivative
    function[out] = ydot(t, y)
        
       out = e*peptide(t)./gsigma(t) - u*y;
        
    end

    % function for g_Sigma levels
    function[out] = gsigma(t)
    
        out = g0 + g1*t;
        
    end


    % function for peptide levels
    function[out] = peptide(t)
    
        out = A*(ktr*t).^(n+1).*exp(-ktr*t)/((n+1)^(n+1.5)*exp(-(n+1)));
        
    end


end