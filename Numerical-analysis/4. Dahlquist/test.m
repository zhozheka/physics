clear;
clc;

N0 = 10;
a = 0;
b = 5;  
global lambda;
lambda = 1;

U0 = 1;
num = 15;

%ERK1
% a11 = 0;

%DIRK1
%a11 = 1;

% a11 = 0.5;

%CROS1
 a11 = (1 + 1i) / 2;


T = zeros(1,N0);

tau = (b-a) / N0;

for i = 1:N0
   T(i) = tau * (i-1); 
end

for k = 1:num
    
    N = N0 * 2^(k-1);
    tau = (b-a) / N;
    u = zeros(N,1);
    t = zeros(N,1);

    t(1) = 0;
    
    
    %первый расчет
   	w = f(U0, tau/2) / (1 - a11*tau*fu(U0, tau/2));
    u(1) = U0 + tau * real(w);
    t(1) = tau;
    
    %считаем на k-й сетке
    for n = 1:N
        
        w = f(u(n), t(n) + tau/2) / (1 - a11*tau*fu(u(n), t(n)+tau/2));
        u(n+1) = u(n) + tau * real(w);
        t(n+1) = t(n) + tau;
    end
    
    %сохраняем в узлах начальной сетки
    for n = 1:N0
        D(k, 1, n) = u((2^(k-1)*n));
    end
    k
    
    for m = 2:k
        p0=2;
        R(k,m-1,:) = ( D(k,m-1,:) - D(k-1,m-1,:) ) / (2^(2+m-2 )  - 1);
        D(k,m,:) = D(k,m-1,:) + R(k,m-1,:);
    end
end

for k = 2:num
   for m = 1:num-1
      P(k,m,:) = log( abs( R(k-1,m,:) ./  R(k,m,:) ) ) / log(2);
   end
end

%в норме лебега

P_norm = P(:,:,1).^2;
for i = 2:N0
    P_norm = P_norm + P(:,:,i).^2;
end
P_norm = sqrt(P_norm / N0);
figure;
hold on;
plot(t,u);
grid on;
hold off;


%{
plot(t,u,'r*-'); 
grid on;
plot(t,q,'b-');
xlabel('t'); 
ylabel('u(t)');
%text(0.2,2.3,'e=0.001, N=6000');
%}
figure;
plot( T, P);
%axis([-0.9 2 1.5 2.5]);
grid on;
