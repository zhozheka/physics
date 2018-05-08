
clear;
clc;

r = [0 0;
    1.52098232*10^11 0;
    1.52098232*10^11+3.84467*10^8 0;
    5.24824*10^12 0];

v = [0 -0.089021687078894;
    0 29270;
    0 29270+1023;
    0 -900];

Z0 = [r, v];

T0 = 365*24*3600*76;
T_ = T0 * 1;

h0 = 0.6e9;

RK = fun(h0, Z0);

n=0;
while ( RK.getT <= T_ )
    n=n+1;
    RK.step(h0);
    Z = RK.getZ;
    T(n) = RK.getT;
    C(n) = RK.getC;
    H(n) = RK.getH;
    tmp = RK.getZ;
    X(n,:) = tmp(:,1);
    Y(n,:) = tmp(:,2);
    Rad(n) = sqrt(X(4)^2 + Y(4)^2);
    RK.getT
end
figure;
hold on;

    plot(X(:,1), Y(:,1), 'go');
    plot(X(:,2), Y(:,2));
    plot(X(:,4), Y(:,4));
    
    legend('Солнце', 'Земля', 'Комета');
    xlabel('X');
    ylabel('Y');
    
    grid on;
hold off;

figure;
hold on;
    plot(T, Rad, 'b-');
    plot(T, C*1e19, 'g-');
    legend('Расстояние от Солнца','Кривизна * 10^{19}')
    xlabel('Время');
    grid on;
hold off;

figure;
hold on;
    plot(T, H);
    xlabel('Время');
    ylabel('Шаг');
    grid on;
hold off;
