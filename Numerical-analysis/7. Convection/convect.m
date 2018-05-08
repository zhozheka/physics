clear;
clc;


M0 = 200;
global N0;
N0 = 50;
a = 0;
b = 1;
global dx;
dx = (b - a) / N0;
x = zeros(1, N0);

for n = 1:N0
    x(n) = a + n * dx;
end

dt = 0.3 / (M0 - 1);
t = zeros(1, M0);


%параметры
%a11 = 0;
%a11 = 1;
%a11 = 1/2;
a11 = (1 + 1i) / 2;

u = zeros(M0, N0);

%начальные условия
for n = 1:N0
    u(1,n) = -x(n) + 1;
end

for j = 1:M0-1
    w = ( (eye(N0) - a11 * dt * fu(u(j,:), t(j)) )^(-1)) * f(u(j,:), t(j)+dt/2);
    u(j+1,:) = u(j,:) + dt * real(w)';
    t(j+1) = t(j) + dt;
end

figure;
surf(x, t, u);

figure;
k=1;
for m=1:k:M0
    plot(x,u(m,:),'-o','MarkerSize',3);
    hold on;
    axis([a b 0 2.5]);
    
    
    hold off;
    drawnow;
    pause(0.0001);
    
end

for j = 1:M0
    P_t(j) =  sqrt(sum(P(j,:).^2) / N0);
end


