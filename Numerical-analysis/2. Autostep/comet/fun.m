 classdef fun < handle   
    properties    
        l, h_, h, R;
        c = zeros(3);
    end
    
    methods
        function this = fun(h0_, Z0)
            
            this.R = [Z0æ 0];
            this.l = 0;
            this.h_  = h0_;
            this.h = h0_;
        end
        
                     
        function step (this, h)
            
            if (h ~= 0)
                this.h = h;
            end         
            
            k1 = this.Rt(this.R);
            k2 = this.Rt(this.R + (this.h/2)*k1 );
            k3 = this.Rt(this.R + (this.h/2)*k2 );
            k4 = this.Rt(this.R + this.h*k3);

            this.R = this.R + ( k1 + k2*2 + k3*2 + k4 ) * (this.h / 6);
            
            this.l = this.l + this.h;
            this.c = (this.Rt(this.R)*3 - k2*2 - k3*2 + k1) / this.h;
            this.h = this.h_ / (1 + 1e4* (this.c*(this.c)')^0.25);
            
        end
        
        
        function[value] = Rt(this, R)
            T = R(5);
            Z = R(1:4);
            f = F_comet(Z);
            Zl= f / sqrt(1 + norm(f));
            Tl = 1 / sqrt(1 + norm(f));
            value = [Zl, Tl];
            
        end
       
        
        function [value] = getZ(this)
            
            value = this.R(1:4); 
        end
        
           
        function [value] = getT(this)
           value = this.R(5); 
        end
        
        function [value] = getL(this)
           value = this.l; 
        end
        
        function [value] = getH(this)
           value = this.h; 
        end
        
        function [value] = getC(this)
           value = norm(this.c);
        end
        
    end
end

