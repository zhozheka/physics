clear;
clc;

%1 алгебр, 12 алг
coord = [10;  0.2; 10; -0.2];
vel =  zeros(4,1);

global C;
global l;
C = 1;
l = sqrt((coord(1) - coord(3))^2 + (coord(2) - coord(4))^2);

u0 =  [coord ; vel; 0]';
u(1,:) = u0;

t(1) = 0;
tau = 0.01;

N = 60. / tau;  

M = eye(9,9);
M(9,9) = 0;
a11 = (1+1i) / 2;
len(1) = l;

for n = 1:N-1

    w = f(u(n,:)) * (M - tau*a11*fu( u(n,:) ) )^(-1);
    u(n+1,:) = u(n,:) + tau * real(w);
    t(n+1) = t(n) + tau;
    len(n+1) = sqrt( (u(n+1,1)-u(n+1,3))^2 + (u(n+1,2)-u(n+1,4))^2);
    
end

figure;
hold on;
title('Диполи');
plot(u(:,1), u(:,2), 'b-' );
plot(u(:,3), u(:,4), 'g-' );
legend('+','-');
xlabel('X');
ylabel('Y');
grid on;
hold off;

figure;
hold on;
title('Расстояние между зарядами');
plot(t, len);
xlabel('время');
ylabel('расстояние');
grid on;
hold off;