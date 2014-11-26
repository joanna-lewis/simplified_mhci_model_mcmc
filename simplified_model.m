function[sol] = simplified_model(npep, parameters, timepoints)

    % simplified model of peptide presentation
    
    % npep is the number of peptides
    
    % parameters contains the model parameters: A, the scaling for each
    % intracellular peptide model; n, the number of transit compartments
    % for each intracellular peptide model; ktr, the transit rate constant
    % between compartments in the intracellular model, and u which captures
    % a combination of MHC-peptide unbinding rate and the scaling between 
    
    % timepoints contain the times of the observations, starting at time 0.
    
    % set up vector for solution
    sol = zeros(2*npep, numel(timepoints));
    
    for i=1:npep
    
        data.p = [parameters(i), ... % A
            parameters(npep + i), ... % n 
            parameters(2*npep + i), ... % ktr
            parameters(3*npep + i), ... % u
            parameters(4*npep + 1) ];
            
        % peptide levels
        sol(i,:) = 0.01*peptide(timepoints);
        
        %%%%%%%%%%%%%%%%%%%%%%%
        % Initialise cvodes
        %%%%%%%%%%%%%%%%%%%%%%%
        
        options = CVodeSetOptions(  'UserData', data,...
                                     'RelTol',1e-10,...
                                     'AbsTol',1e-5,...
                                    ...%'MaxStep',600, ...
                                    'LinearSolver','Dense',...
                                    'MaxNumSteps', 5000);
                
        CVodeInit(@ydot, 'BDF', 'Newton', 0, 0, options);

        %%%%%%%%%%%%%%%%%%%%%%%
        % solve odes
        %%%%%%%%%%%%%%%%%%%%%%%

        sol(npep+i,1) = 0;
        [~, ~, sol(npep+i,2:end)] = CVode(timepoints(2:end), 'Normal');

        CVodeFree;
        
    end
    
    % function for calculating derivative
    function[y_dot, flag, new_data] = ydot(t, y, data)

        y_dot = 1e4*peptide(t)./gsigma(t) - data.p(4)*y;

        flag=0;
        new_data = [];

    end

    % function for g_Sigma levels
    function[out] = gsigma(t)
    
        out = 1 + data.p(5)*t;
        
    end


    % function for peptide levels
    function[out] = peptide(t)
    
        out = data.p(1)*(data.p(3)*t).^(data.p(2)+1).*exp(-data.p(3)*t)/((data.p(2)+1)^(data.p(2)+1.5)*exp(-(data.p(2)+1)));
        
    end


end