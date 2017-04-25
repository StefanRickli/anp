classdef polyratio < handle
    properties
        num
        denom
    end
    
    methods
        function this = polyratio(a,b)
            assert(isa(a,'double') && isa(b,'double'));
            
            this.num = a;
            this.denom = b;
        end
        
        function y = normalize(this)
            % Pads the shorter of num and denum with zeros such that they are of equal length.
            
            [this.num,this.denom] = pad_left(this.num,this.denom);
            y = this;
        end
        
        function varargout = reduce(this)
            % ATTENTION: This method is dangerous as it deletes poles and
            %            zeros that are close enough ( < tol) to each other
            
            % TODO: tol most certainly needs to be dynamically changed
            %       based on the expected precision of 'roots'
            tol = 1e-4;
            
            zeros = roots(this.num);
            poles = roots(this.denom);
            
            while(true)
                dist = inf;
                zi = NaN;
                pi = NaN;
                
                for ii = 1:length(zeros)
                    for jj = 1:length(poles)
                        e = abs(zeros(ii) - poles(jj));
                        if e < dist
                            dist = e;
                            zi = ii;
                            pi = jj;
                        end
                    end
                end
                
                if (dist < tol)
                    zeros(zi) = [];
                    poles(pi) = [];
                else
                    break;
                end
            end
            
            this.num =      poly(zeros);
            this.denom =    poly(poles);
            
            this.normalize;
            
            switch nargout
                case 0
                    varargout = {};
                case 1
                    varargout{1} = this;
                case 2
                    varargout{1} = zeros;
                    varargout{2} = poles;
                otherwise
                    error('''reduce'' accepts either one output argument, returning the object itself, or two output arguments, returning zeros and poles of the polyratio');
            end
        end
        
        function y = mult(this,x)
            assert(isa(x,'polyratio'));
            
            y = polyratio(conv(this.num,x.num),conv(this.denom,x.denom));
        end
        
        function y = add(this,x)
            assert(isa(x,'polyratio'));
            
            [a1,a2] = pad_left(conv(this.num,x.denom),conv(x.num,this.denom));
            
            common_denom = conv(this.denom,x.denom);
            
            y = polyratio(a1 + a2, common_denom);
        end
        
        function y = sub(this,x)
            assert(isa(x,'polyratio'));
            
            [a1,a2] = pad_left(conv(this.num,x.denom),conv(x.num,this.denom));
            
            common_denom = conv(this.denom,x.denom);
            
            y = polyratio(a1 - a2, common_denom);
        end
    end
end