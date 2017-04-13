function [x_span,y_span] = get_value_span(~, values) % ignored parameter is 'this'
    
    real_values = real(values);
    imag_values = imag(values);
    
    x_span = [min(real_values),max(real_values)];
    y_span = [min(imag_values),max(imag_values)];


    width = diff(x_span);
    height = diff(y_span);

    % make this plot also rectangular
    if width >= height
        y0 = mean(y_span);
        y_span = [(y0-width/2),(y0+width/2)];
    else
        x0 = mean(x_span);
        x_span = [(x0-height/2),(x0+height/2)];
    end
end