clear;
clc;

global N h;

a11 = (1 + 1i) / 2;

N0 = 20; %��������� ����� ����� �� x
M0 = 50; %��������� ����� ����� �� t
T0 = 0.3; %����� �������


%�� ���������� ����� (��������� � ���������)
tau = T0 / (M0);
h = 1.0 / N0;
    
for n = 1:N0
    x0(n) = n*h;
end

for j = 1:M0
    t0(j)  = j*tau;
end

%������� 3 ����
for k = 1:3
    clear u;
    M = M0 * 2^(k-1); %� 2 ���� �� �������
    N = N0 * 4^(k-1); %� 4 ���� �� ����������
    
    tau = T0 / (M);
    h = 1.0 / N;
    
    %����� �� x
    for n = 1:N
        x(n) = n*h;
    end

    %��������� �������
    for n = 1:N
        u0(n) = -x(n) + 1;
    end
    
    %������ ��� � �������� ����
    w = ( (eye(N) - a11 * tau * fu(u0, 0) )^(-1)) * f(u0, 0+tau/2);
    u(1,:) = u0 + tau * real(w)';
    t(1) = tau;
    
    
    %� �� ���������
    for j = 1:M-1
        w = ( (eye(N) - a11 * tau * fu(u(j,:), t(j)) )^(-1)) * f(u(j,:), t(j)+tau/2);
        u(j+1,:) = u(j,:) + tau * real(w)';
        t(j+1) = t(j) + tau;
    end
    
    
    %������ ������� �� ������� �����
    figure;
    hold on;
        surf(x, t, u, 'EdgeColor','none','LineStyle','none','FaceLighting','phong');
        grid on;
        title (['������� ��� N = ' num2str(N) ', M = ', num2str(M)]);
        xlabel('x');
        ylabel('t');
    hold off;
    
    
    %��������� �� ��������� �����
    for j = 1:M0
        for i = 1:N0
            D(k, j, i) = u(j*2^(k-1), i*4^(k-1));
        end
    end
    
end   
   

%������� ����������� ������� �������� � ������ �����
P(:,:) = log( abs(( D(2, :, :) - D(1, :, :)) ./ (D(3, :, :) - D(2, :, :)) ) ) / log(4);

%������ ����������� ������� �� �������
figure;
k=1;
for m=1:k:M0
    plot(x0,P(m,:),'-o','MarkerSize',3);
    hold on;
        axis([0 1 -1 3]);
        title('������� �������� � ������� �������');
        xlabel('x');
        ylabel('������� ��������');
        grid on;
        
    hold off;
    drawnow;
    pause(0.01);  
end

% � ��������� ������ ������� � ����� �=1 ������� ������� �� 1, ��-�� ����, 
% ��� �� ������ ������� ������� �������� �����.


for j = 1:M0
    P_t(j) =  sqrt(sum(( P(j,:).^2)) / N0);
end
figure;
plot(t0, P_t)

