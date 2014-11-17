function out = LL(pred, tryptic_data, Kb_data, error)

%error = 0.06;

    %error = [0.01, 0.02, 5, 50];
%     error = [0.02, 0.02, 0.02,  0.2, ... % tryptic error
%                 50, 200, 500, 50 ... % epitope erroe
%                 ];


    [npep, N] = size(tryptic_data);
    
    %error(2) = 50;

    out = 0;
    
    for i=1:npep
    
        out = out ...
                - 0.5*N*log(2*pi*error(i)^2)  ...
                - 0.5*(tryptic_data(i,:) - pred(i,:))*(tryptic_data(i,:) - pred(i,:))'/error(i)^2 ...
                - 0.5*N*log(2*pi*error(npep+i)^2)  ...
                - 0.5*(Kb_data(i,:) - pred(npep+i,:))*(Kb_data(i,:) - pred(npep+i,:))'/error(npep+i)^2;
    end
            
end