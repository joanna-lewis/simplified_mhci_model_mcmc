function out = LL(pred, tryptic_data, Kb_data, error)

%error = 0.06;

    N = numel(pred);
    
    %error(2) = 50;

    out = - 0.5*N*log(2*pi*error(1)^2)  ...
            - 0.5*(tryptic_data - pred(1,:))*(tryptic_data - pred(1,:))'/error(1)^2 ...
            - 0.5*N*log(2*pi*error(2)^2)  ...
            - 0.5*(Kb_data - pred(2,:))*(Kb_data - pred(2,:))'/error(2)^2;

end