clear;
clc;

global a_ l0_ h0 u0 w_ sigma_ ;


for k = k_start : k_max
    
    dl = h0 / 2^(k-1);
    RK = fun(dl);
    delta = 0;
    
    %нулевой шаг
    n = 0;
    H0(1)=RK.getH;
    U0(1)=RK.getU;
    T0(1)=RK.getT;
    
    while ( RK.getL < L0 )
        
        n = n + 1;
        RK.step(0); %делаем шаг
        H(1,n) = RK.getH;     
        L(k,n) = RK.getL;
        U(1,n) = RK.getU;
        T(1,n) = RK.getT;
        
        %начиная со 2го, только для четных
        if (k > k_start) && (rem(n,2)==0) && (n <= N(k-1))
            delta = delta + (L(k,n) - L(k-1,n/2))^2;            
        end
    end
    k
    
    Delta(k) = sqrt( 2*delta / n);
    N(k) = n;
%     DELTA = 0.003;
    if (Delta(k) < DELTA)  && (k > 1)
        break;
    end
end


% figure;
%     hold on
%     plot(Delta, 'bo-');
%     title(['Сравниваем с критерием \delta = ', num2str(DELTA), ', \lambda = ' num2str(l0_)]);
%     xlabel('Номер сгущения');
%     ylabel('Сравниваем с \delta')
%     
% grid on;

%второй этап-----------------------

ERR_real = zeros(2,1);
ERR_est = zeros(2,1);
ERR_real_C = zeros(2,1);
ERR_est_C = zeros(2,1);
N = zeros(2,1);

% num = 8;
N(1) = n;
L2(1,:) =  L(k,:);

for k = 2:num
  
    %пересчитываем шаги на k-ом сгущении
    %первые шаги
    H0(k)  = H0(k-1) * sqrt(H0(k-1))  / (sqrt(H0(k-1)) + sqrt(H(k-1,1)));
    H(k,1) = H0(k-1) * sqrt(H(k-1,1)) / (sqrt(H0(k-1)) + sqrt(H(k-1,1)));
    
    H(k,2) = H(k-1,1) * H0(k-1)^0.25  / (H0(k-1)^0.25 + H(k-1, 2)^0.25);
    H(k,3) = H(k-1,1) * H(k-1,2)^0.25 / (H0(k-1)^0.25 + H(k-1, 2)^0.25);
    
    for s = 2:n-2
        
        H(k,2*s)  =  H(k-1,s) * H(k-1,s-1)^0.25 / (H(k-1,s+1)^0.25 + H(k-1,s-1)^0.25);
        H(k,2*s+1) = H(k-1,s) * H(k-1,s+1)^0.25 / (H(k-1,s+1)^0.25 + H(k-1,s-1)^0.25);       
    
    end
    %два последних шага
    H(k,2*n-2) = H(k-1,n-1) * sqrt(H(k-1,n-2)) / (sqrt(H(k-1,n-2)) + sqrt(H(k-1,n-1)));
    H(k,2*n-1) = H(k-1,n-1) * sqrt(H(k-1,n-1)) / (sqrt(H(k-1,n-2)) + sqrt(H(k-1,n-1)));
    
    %возвращаемся в начало    
    RK2 = fun(0);
    
    %первая точка
    U0(k) = RK2.getU;
    T0(k) = RK2.getT;
    
    %первый шаг
    RK2.step(H0(k));
    U(k,1) = RK2.getU;
    T(k,1) = RK2.getT;
    L2(k,1) = RK2.getL;
    C(1) = RK2.getCappa;
    
    N(k) = N(k-1)*2;
    
    n = 2*n; %число шагов на текущем сгущении
    for s = 2:n
        RK2.step(H(k,s-1));
        U(k,s) = RK2.getU;
        T(k,s) = RK2.getT;
        L2(k,s) = RK2.getL;
        C(s) = RK2.getCappa;
    end
    k
    %теперь считаем погрешности:----------------------------
    %РЕАЛЬНАЯ       
    
    ERR_real(k) = 0;
    ERR_real_C(k) = 0;
    e1=zeros(2,1);
    j=0;
    for s = 1:n
        j=j+1;
        ERR_real(k) = ERR_real(k) + (U(k,s) - U_real(T(k,s)))^2;
        
        e1(s) = abs(U(k,s) - U_real(T(k,s)));
        if e1(s) > ERR_real_C(k)
             ERR_real_C(k) = e1(s);
        end
    end
    
    ERR_real(k) = sqrt(ERR_real(k)/j);
    
    
    %ОЦЕНКАdelta
    
    ERR_est(k) = 0;
    ERR_est_C(k) = 0;
    e2=zeros(2,1);

    for s = 1:N(k-1)
        
        eT = (T(k-1,s) - T(k,s*2) ) / 15;
        eU = (U(k-1,s) - U(k,s*2) ) / 15;

        ERR_est(k) = ERR_est(k) + (eU - eT*Ut( T(k,s*2), U(k,s*2)))^2;
        
        e2(s) = abs(eU - eT*Ut( T(k,s*2), U(k,s*2)));

        if e2(s) > ERR_est_C(k)
            ERR_est_C(k) = e2(s);
        end
 
    end   
    ERR_est(k) = sqrt( ERR_est(k) / (N(k-1)-2) );
end

%графики погрешностей
% figure;
    hold on;    
    title (['Погрешность для задачи с \lambda = ' num2str(l0_), ', \delta = ' num2str(DELTA) ]);
    

    plot ( log10(N), log10(ERR_est), '-b*');
    plot ( log10(N), log10(ERR_real), '-g*');
%     
%     plot ( log10(N), log10(ERR_est_C), '-r*');
%     plot ( log10(N), log10(ERR_real_C), '-k*');
%     
    legend('Оценка ошибки', 'Реальная ошибка', 'Оценка ошибки (C)', 'Реальная ошибка (C)');

    xlabel ('lg(N)');
    ylabel ('lg(Error)');
    grid on;
    hold off;

%структура решения
% figure;
% hold on;
% title (['Структура решения для задачи с \lambda = ' num2str(l0_)]);
% plot(T(k,:), U(k,:));
% plot(T(k,:), C(:));
% xlabel ('T');
% ylabel ('U, \kappa');
% legend('Вид решения', 'Кривизна');
% 
% grid on;
% hold off;
