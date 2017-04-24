function [c,d] = pad_left(a,b)
    % Pads the shorter vector on the left with zeros such that both vectors have the same length.
        dl = length(a) - length(b);
        
        c = [zeros(1,max(0,-dl)),a];
        d = [zeros(1,max(0, dl)),b];
end
