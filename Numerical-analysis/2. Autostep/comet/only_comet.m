
clear;
clc;
r = [5.24824*10^12 0]; 
v = [0 -900];
Z0 = [r, v];

T0 = 365*24*3600*76;
T_ = T0 * 1;

h0 = 1e9;

RK = fun(h0, Z0);

n=0;
while ( RK.getT <= T_ )
    n=n+1;
    RK.step(0);
    Z = RK.getZ;
    X(n) = Z(1);
    Y(n) = Z(2);
    Rad(n) = sqrt( X(n)^2 + Y(n)^2);
    T(n) = RK.getT;
    C(n) = RK.getC;
    Vx(n) = Z(3);
    Vy(n) = Z(4);
    V(n) = sqrt( Vx(n)^2 + Vy(n)^2); 
    H(n) = RK.getH;
    RK.getT
end
figure;
hold on;
    plot(X, Y, 'bo', 'MarkerSize',3);
    plot(0, 0, 'rp', 'MarkerFaceColor', 'y', 'MarkerSize',20);
    legend('Земля', 'Солнце');
    xlabel('X');
    ylabel('Y');
    
    grid on;
hold off;

figure;
hold on;
    plot(T, Rad, 'b-');
    plot(T, C*1e19, 'g-');
    plot(T, V*1e8, 'r-');
    plot(T, H*1e4, 'm-');
    legend('Расстояние от Солнца','Кривизна * 10^{19}', 'Скорость * 10^8', 'Шаг')
    xlabel('Время');
    grid on;
hold off;


