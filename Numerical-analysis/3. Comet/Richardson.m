clear;
clc;

 num = 4;
% m = [1.9891*10^30; 
%     5.9736*10^24;
%     7.3477*10^22; 
%     2.2*10^14]; 

m = [1.9891*10^30; 
     5.9736*10^24;
     7.3477*10^22;
     2.2*10^14]; 

r = [0 0;
    1.52098232*10^11 0;
    1.52098232*10^11+3.84467*10^8 0;
    5.24824*10^12 0]; 
v = [0 -0.089021687078894;
    0 29270;
    0 29270+1023;
    0 -900];
 
Z = [r, v];

X1 = zeros(1000,1);   
Y1 = zeros(1000,1);
N0 = 2000;
%dt0 = 3600*24*10;
T0 = 3600*24*365*75.3;
T= T0 *3
n=1;

%            -----Matrix  Z-----
%
%            |Coords|  |Velocities|
%   | Sun |    *   *      *    *
%   |Earth|    *   *      *    *
%   | Moon|    *   *      *    *
%   |Comet|    *   *      *    *

for k=1:num
    N = N0 * 2^(k-1); %����� �����
    dt = T / N; %��� �� �������
    
    Z = [r, v]; %��������� �������
    
    for n = 1:N
        
        k1 = F(Z);
        k2 = F(Z + k1 * (dt/2) );
        k3 = F(Z + k2 * (dt/2) );
        k4 = F(Z + k3 * dt);
        Z = Z + ( k1 + k2*2 + k3*2 + k4 ) * (dt / 6);        
        
        C(n,:,:) = Z; %��������� �������� ������� �� ������ ��������
    end
    
    %��������� �� ���������� �����
    %D(k,n,:,:) k-����� ��������, n-����� ���� �� ������ �����
    
    n
    for n = 1:N0
        j = n * 2^(k-1);
        D(k, 1, n, :, :) = C(n * 2^(k-1), :, :);
    end
    
    hs(k) = subplot(2,2,k);
        hold on;
        grid on;
        %title([ '����� ����� N = ', num2str(N), ' �������� #', num2str(k)]);
        title ([ '#', num2str(k), ', N = ', num2str(N)]);
        xlabel('X');
        ylabel('Y');
        plot( C(:, 1,1), C(:, 1,2));
        plot( C(:, 2,1), C(:, 2,2)); 
        plot( C(:, 3,1), C(:, 3,2));
        plot( C(:, 4,1), C(:, 4,2)); 
        legend('������', '�����','����','������ ������');
    
    p=4+num;
    
    X(k,1) = Z(4,2);        
    
    for m = 2:k
        %R(k,m-1) = ( X(k,m-1) - X(k-1,m-1)) / (2^p - 1); 
        %X(k,m) = X(k,m-1) + R(k,m-1);
        p1= p + (m-1);
        R(k,m-1,:,:,:) = ( D(k, m-1, :, :, :) - D(k-1, m-1, :, :, :)) / (2^(4+num-2) -1);
        D(k,m,:,:,:) = D(k,m-1,:,:,:) + R(k,m-1,:,:,:);
    end
end

for k = 2:num
    for m = 1:num-1
        %P(k,m) = log(abs(R(k-1,m))/abs(R(k,m)))/log(2);
        P(k,m,:,:,:) = log( abs( R(k-1, m,:,:,:) ./ R(k,m,:,:,:) ) ) / log(2);
    end
end

figure;
hold on;
grid on;
title(['���������']);
C1(:,:,:) = D(3,3,:,:,:);

plot( C1(:, 2,1), C1(:, 2,2)); 
plot( C1(:, 4,1), C1(:, 4,2)); 

legend('�����','������ ������');