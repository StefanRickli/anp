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
            
            last_nonzero_idx = find(this.num,1,'last');
            n_roots_at_zero = length(this.num) - last_nonzero_idx;
            z = zeros(n_roots_at_zero,1);
            zeros_mult = multroot(this.num(1:end-n_roots_at_zero));
            
            last_nonzero_idx = find(this.denom,1,'last');
            n_roots_at_zero = length(this.denom) - last_nonzero_idx;
            p = zeros(n_roots_at_zero,1);
            poles_mult = multroot(this.denom(1:end-n_roots_at_zero));
            
            if ~isempty(zeros_mult)
                for ii = 1:length(zeros_mult(:,1))
                    z = [z;ones(zeros_mult(ii,2),1) * zeros_mult(ii,1)];
                end
            end
            if ~isempty(poles_mult)
                for ii = 1:length(poles_mult(:,1))
                    p = [p;ones(poles_mult(ii,2),1) * poles_mult(ii,1)];
                end
            end
            
            while(true)
                dist = inf;
                zii = NaN;
                pii = NaN;
                
                % Find smallest pairwise distance between poles and zeros
                for ii = 1:length(z)
                    for jj = 1:length(p)
                        e = abs(z(ii) - p(jj));
                        if e < dist
                            dist = e;
                            zii = ii;
                            pii = jj;
                        end
                    end
                end
                
                if (dist < tol)
                    fprintf('Removing zero (%s) and pole (%s) from polynomial fraction. Distance is: %s\n', num2str(z(zii)), num2str(p(pii)), num2str(abs(z(zii) - p(pii))));
                    z(zii) = [];
                    p(pii) = [];
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
            this.num =      real(num_scale * poly(z));
            this.denom =    real(denom_scale * poly(p));
            
            switch nargout
                case 0
                    varargout = {};
                case 1
                    varargout{1} = this;
                case 2
                    varargout{1} = z;
                    varargout{2} = p;
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
