function dudy = cheb_eval_series(b, x)
    b = b(:);
    x = x(:);
    N = length(b)-1;

    T0 = ones(size(x));
    if N==0
        dudy = b(1)*T0;
        return
    end

    T1 = x;
    dudy = b(1)*T0 + b(2)*T1;

    for n = 1:N-1
        T2 = 2*x.*T1 - T0;       % T_(n+1)
        dudy = dudy + b(n+2)*T2;
        T0 = T1;
        T1 = T2;
    end
end