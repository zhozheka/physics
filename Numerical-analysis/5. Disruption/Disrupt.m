clear;
clc;

N0 = 1000;
a = 0;
b = 5;  
global lambda;
lambda = 1;

U0 = 1;
num = 5;

%ERK1
%a11 = 0;

%DIRK1
%a11 = 1;

%a11 = 0.5;

%CROS1
a11 = (1 + 1i) / 2;


U = zeros(num,N0);
R = zeros(num,N0);
P = zeros(num,N0);
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

    u(1) = 1;
    t(1) = 0;
    
    %считаем на k-й сетке
    for n = 1:N-1
        w = f(u(n), t(n) + tau/2) / (1 - a11*tau*fu(u(n), t(n)+tau/2));
        u(n+1) = u(n) + tau * real(w);
        t(n+1) = t(n) + tau;
    end
    
    %сохраняем в узлах начальной сетки
    for n = 1:N0
        U(k, n) = u((2^(k-1)*(n-1)+1));
    end
    
end

figure;
hold on;
plot(t,u);
grid on;
hold off;

for k = 2:num
    for n=1:N0
        R(k,n) = ( U(k,n) - U(k-1,n) ) / 3;
    end
end

for k = 3:num
   for n = 2:N0
       P(k,n) = abs(log(R(k-1,n)/R(k,n)))/log(2);
   end
end
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
