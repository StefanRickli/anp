classdef tools
    methods(Static)
        function y = map(x,t0,t1,u0,u1)
            y = ((u0-u1).*x + (t0*u1-t1*u0))/(t0-t1);
        end
        
        function y = iterator_modulo(x,m)
            y = mod(x - 1,m) + 1;
        end
        
        function dbg(varargin)
            global debug_text
            if debug_text
                fprintf(varargin{:});
            end
        end
    end
end