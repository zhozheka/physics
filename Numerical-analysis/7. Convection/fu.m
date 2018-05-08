function fu = fu(U, t)
    global N h;
    fu = zeros(N, N);

        fu(1, 1) = -(2*U(1) - exp(-t))/h + 2*U(1)*exp(U(1)^2);
    for n = 2:N
        fu(n, n) = -(2*U(n) - U(n-1))/h + 2*U(n)*exp(U(n)^2);
        fu(n, n - 1) = U(n)/h;
    end
end

