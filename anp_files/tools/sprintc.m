function str = sprintc(z)
    % Converts a complex number into a string
    y = imag(z);
    if y >= 0
        str = sprintf('%.1f+%.1fi', real(z), y);
    else
        str = sprintf('%.1f-%.1fi', real(z), -y);
    end
end