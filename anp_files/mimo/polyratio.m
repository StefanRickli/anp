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
            
            this.normalize();
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
            
            fprintf('polyratio.reduce: Will remove all pole/zero combinations that have a difference of < %s\n', num2str(tol,'%1.1e'));
            fprintf('This is a (quick and) dirty numerical solution. Keep this in mind!\n');
            
            % Handle constant polynomials and the zero function outside uvFactor
            if isempty(find(this.num ~= 0, 1, 'first')) || find(this.num ~= 0, 1, 'first') == length(this.num)
                zeros = [];
            else
                % uvFactor doesn't like vectors with leading zeros. Get rid
                % of them.
                zeros_res =	uvFactor(this.num(find(this.num ~= 0, 1, 'first'):end));
                zeros = uvFactor_res2roots(zeros_res);
            end
            
            if isempty(find(this.denom ~= 0, 1, 'first')) || find(this.denom ~= 0, 1, 'first') == length(this.denom)
                poles = [];
            else
                poles_res =	uvFactor(this.denom(find(this.denom ~= 0, 1, 'first'):end));
                poles = uvFactor_res2roots(poles_res);
            end
            
            while(true)
                dist = inf;
                zii = NaN;
                pii = NaN;
                
                % Find smallest pairwise distance between poles and zeros
                for ii = 1:length(zeros)
                    for jj = 1:length(poles)
                        e = abs(zeros(ii) - poles(jj));
                        if e < dist
                            dist = e;
                            zii = ii;
                            pii = jj;
                        end
                    end
                end
                
                if (dist < tol)
                    fprintf('Removing zero (%s) and pole (%s) from polynomial fraction. Distance is: %s\n', num2str(zeros(zii)), num2str(poles(pii)), num2str(abs(zeros(zii) - poles(pii))));
                    zeros(zii) = [];
                    poles(pii) = [];
                else
                    % The smallest distance between any pole/zero
                    % combination is higher than our threshold, 'tol'. This
                    % means that we're done with polynomial reduction.
                    break;
                end
            end
            
            % We need to know the coefficient of the highest power of the
            % polynomials apart from their roots, as we loose this constant
            % otherwise.
            num_scale =      this.num(find(this.num ~= 0,1,'first'));
            denom_scale =    this.denom(find(this.denom ~= 0,1,'first'));
            
            % Put the numerator and denominator together again.
            this.num =      num_scale * poly(zeros);
            this.denom =    denom_scale * poly(poles);
            
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
            
            y.reduce;
            y.normalize;
        end
        
        function y = add(this,x)
            assert(isa(x,'polyratio'));
            
            [a1,a2] = pad_left(conv(this.num,x.denom),conv(x.num,this.denom));
            
            common_denom = conv(this.denom,x.denom);
            
            y = polyratio(a1 + a2, common_denom);
            
            y.reduce;
            y.normalize;
        end
        
        function y = sub(this,x)
            assert(isa(x,'polyratio'));
            
            [a1,a2] = pad_left(conv(this.num,x.denom),conv(x.num,this.denom));
            
            common_denom = conv(this.denom,x.denom);
            
            y = polyratio(a1 - a2, common_denom);
            
            y.reduce;
            y.normalize;
        end
    end
end

function poly_roots = uvFactor_res2roots(factors)
    % This function converts the output of uvFactor into a column vector of roots
    
    poly_roots = [];
    [m,~] = size(factors);
    
    for ii = 1:m
        for jj = 1:factors(ii, 3)
            poly_roots(end+1) = -factors(ii,2)./factors(ii,1);
        end
    end
    
    poly_roots = poly_roots';
end