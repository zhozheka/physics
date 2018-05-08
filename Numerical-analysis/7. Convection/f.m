function f = f(U, t)
    global N h;
    f = zeros(1,N);

    f(1) = - U(1)*(U(1) - exp(-t)) / h + exp(U(1)^2);
    for n = 2: N
        f(n) = - U(n)*(U(n) - U(n-1)) / h + exp(U(n)^2);
    end
    f = f';

end

