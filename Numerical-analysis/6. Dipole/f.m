function [ f ] = f( u )
    %F הטפפ קאסעü
    global C l;
    
    x1 = u(1);
    y1 = u(2);
    x2 = u(3);
    y2 = u(4);
    
    v1x = u(5);
    v1y = u(6);
    v2x = u(7);
    v2y = u(8);
    
    
    f(1) = v1x; %v1_x
    f(2) = v1y; %v1_y
    f(3) = v2x; %v2_x
    f(4) = v2y; %v2_y
    
    r1 = sqrt(x1^2 + y1^2);
    r2 = sqrt(x2^2 + y2^2);
    
    a1x = C*x1 / r1^3;
    a1y = C*y1 / r1^3;
    a2x = -C*x2 / r2^3;
    a2y = -C*y2 / r2^3;
    
    
    f(5) = a1x;
    f(6) = a1y;
    f(7) = a2x;
    f(8) = a2y;
    
    f(9) = (x1 - x2)^2 + (y1 - y2)^2 - l^2;
    
end

