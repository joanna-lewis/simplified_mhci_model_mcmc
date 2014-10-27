function[sol] = simplified_model_tryptic(parameters, timepoints)

    % simplified model of peptide presentation
    
    % allocate parameters
    A = parameters(1);
    n = parameters(2);
    ktr = parameters(3);
    
%     e = parameters(1);
%     
%     g0 = parameters(1);
%     g1 = parameters(1);
%     
%     u = parameters.u;
%     
    % ode solver options
    options = odeset('RelTol', 10^-10, 'AbsTol', 10^-10);
    
    % set up vector for solution
    sol = zeros(1, numel(timepoints));
    
    % peptide levels
    sol(1,:) = peptide(timepoints);
    
    %[~, sol(2,:)] = ode15s(@ydot, timepoints, 0, options);

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
    
        out = A*ktr*(ktr*t).^n.*exp(-ktr*t)/(n^(n+0.5)*exp(-n));
        
    end


end