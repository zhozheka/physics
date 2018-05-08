clear;
clc;

global N h;

a11 = (1 + 1i) / 2;

N0 = 20; %начальное число узлов по x
M0 = 50; %начальное число узлов по t
T0 = 0.3; %время расчета


%на глобальной сетке (совпадает с начальной)
tau = T0 / (M0);
h = 1.0 / N0;
    
for n = 1:N0
    x0(n) = n*h;
end

for j = 1:M0
    t0(j)  = j*tau;
end

%сгущаем 3 раза
for k = 1:3
    clear u;
    M = M0 * 2^(k-1); %в 2 раза по времени
    N = N0 * 4^(k-1); %в 4 раза по координате
    
    tau = T0 / (M);
    h = 1.0 / N;
    
    %сетка по x
    for n = 1:N
        x(n) = n*h;
    end

    %начальные условия
    for n = 1:N
        u0(n) = -x(n) + 1;
    end
    
    %делаем шаг с нулевого слоя
    w = ( (eye(N) - a11 * tau * fu(u0, 0) )^(-1)) * f(u0, 0+tau/2);
    u(1,:) = u0 + tau * real(w)';
    t(1) = tau;
    
    
    %и на остальных
    for j = 1:M-1
        w = ( (eye(N) - a11 * tau * fu(u(j,:), t(j)) )^(-1)) * f(u(j,:), t(j)+tau/2);
        u(j+1,:) = u(j,:) + tau * real(w)';
        t(j+1) = t(j) + tau;
    end
    
    
    %строим решение на текущей сетке
    figure;
    hold on;
        surf(x, t, u, 'EdgeColor','none','LineStyle','none','FaceLighting','phong');
        grid on;
        title (['Решение для N = ' num2str(N) ', M = ', num2str(M)]);
        xlabel('x');
        ylabel('t');
    hold off;
    
    
    %сохраняем на начальную сетку
    for j = 1:M0
        for i = 1:N0
            D(k, j, i) = u(j*2^(k-1), i*4^(k-1));
        end
    end
    
end   
   

%считаем эффективный порядок точности в каждой точке
P(:,:) = log( abs(( D(2, :, :) - D(1, :, :)) ./ (D(3, :, :) - D(2, :, :)) ) ) / log(4);

%строим зависимость порядка от времени
figure;
k=1;
for m=1:k:M0
    plot(x0,P(m,:),'-o','MarkerSize',3);
    hold on;
        axis([0 1 -1 3]);
        title('Порядки точности в течение времени');
        xlabel('x');
        ylabel('Порядок точности');
        grid on;
        
    hold off;
    drawnow;
    pause(0.01);  
end

% В начальный момент времени в около х=1 порядки отличны от 1, из-за того, 
% что на правой границе функция известна точно.


for j = 1:M0
    P_t(j) =  sqrt(sum(( P(j,:).^2)) / N0);
end
figure;
plot(t0, P_t)

