function out = LL(pred, tryptic_data, Kb_data, error)

    % pred are simulations from current parameter values
    % tryptic_data is intracellular peptide levels
    % Kb_data is cell-surface epitope levels
    % error is a two-element vector, [error in tryptic level, error in
    % epitope level] Both elements are standard deviation of measurement
    % error

    [npep, N] = size(tryptic_data);
    
    out = 0;
    
    for i=1:npep
    
        out = out ...
                - 0.5*N*log(2*pi*error(i)^2)  ...
                - 0.5*(tryptic_data(i,:) - pred(i,:))*(tryptic_data(i,:) - pred(i,:))'/error(i)^2 ...
                - 0.5*N*log(2*pi*error(npep+i)^2)  ...
                - 0.5*(Kb_data(i,:) - pred(npep+i,:))*(Kb_data(i,:) - pred(npep+i,:))'/error(npep+i)^2;
    end
            
end