classdef polyratio < handle
    properties
        num
        denom
    end
    
    methods
        function this = polyratio(a,b)
            assert(isequal(class(a),'double') && isequal(class(b),'double'));
            
            this.num = a;
            this.denom = b;
        end
        
        function y = normalize(this)
            % Pads the shorter of num and denum with zeros such that they are of equal length.
            
            [this.num,this.denom] = pad_left(this.num,this.denom);
            y = this;
        end
        
        function y = mult(this,x)
            assert(isequal(class(x),'polyratio'));
            
            y = polyratio(conv(this.num,x.num),conv(this.denom,x.denom));
        end
        
        function y = add(this,x)
            assert(isequal(class(x),'polyratio'));
            
            [a1,a2] = pad_left(conv(this.num,x.denom),conv(x.num,this.denom));
            
            common_denom = conv(this.denom,x.denom);
            
            y = polyratio(a1 + a2, common_denom);
        end
        
        function y = sub(this,x)
            assert(isequal(class(x),'polyratio'));
            
            [a1,a2] = pad_left(conv(this.num,x.denom),conv(x.num,this.denom));
            
            common_denom = conv(this.denom,x.denom);
            
            y = polyratio(a1 - a2, common_denom);
        end
    end
end