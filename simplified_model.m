function[sol] = simplified_model(parameters, timepoints)

    % simplified model of peptide presentation
    
    % allocate parameters
    A = parameters(1);
    n = parameters(2);
    ktr = parameters(3);
    
    e = 1e4;
    
    g0 = 1;
    g1 = parameters(4); %6.74;
    
    u = parameters(5); %5.9327e+01; %3600*10^-4;
         
    % ode solver options
    options = odeset('RelTol', 10^-10, 'AbsTol', 10^-10);
    
    % set up vector for solution
    sol = zeros(1, numel(timepoints));
    
    % peptide levels
    sol(1,:) = 0.01*peptide(timepoints);
        
    [~, sol(2,:)] = ode15s(@ydot, timepoints, 0, options);
    sol(2,:) = sol(2,:);
    
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