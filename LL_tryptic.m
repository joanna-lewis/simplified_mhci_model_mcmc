function out = LL_tryptic(pred, data, error)

    N = numel(pred);

    out = -0.5*N*log(2*pi*error^2) - 0.5*(data - pred)*(data - pred)'/error^2;

end