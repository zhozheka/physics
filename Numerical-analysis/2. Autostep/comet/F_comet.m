function [ Zt ] = F_comet( Z )
    %R: x, y, Vx, Vy
    
    M = 1.9891*10^30;
    G = 6.67300*10^(-11);
    
    x = Z(1);
    y = Z(2);
    r = sqrt(x^2 + y^2);
    
    %Rt = zeros(1,4);
    Zt(1) = Z(3);
    Zt(2) = Z(4);
    Zt(3) =  -G*M*x / r^3;
    Zt(4) =  -G*M*y / r^3;
end

