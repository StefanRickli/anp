function z = test_evaluate_input_function(t,poles,zeros,radii,angles)
    if any(strcmp(who('global'),'debug'))
        global debug;
    else
        debug = false;
    end
    
%     poles = 1i*[-2,-1,-0.6,0.6,1,2];
%     zeros = 1i*[0];
%     radii.inf = 10;
%     angles.crop = 7*pi/180;
%     angles.detour = 45*pi/180;    
    separation_max_pole = 1/4;
    separation_max_zero = 1/8;

    interval_list = test_fill_interval_list_t_to_q(poles,zeros,radii,angles,separation_max_pole,separation_max_zero);
    
%     q_max = interval_list(end).q(2);
%     t = 0:1/1000:1;
%     q = map(t,0,1,0,q_max);
%     z = evaluate_interval_q(interval_list,q);
%     t = 0:1/1000:1;
    
    z = evaluate_interval_t(interval_list,t(1:end));
    
    if debug
        figure;
        scatter(real(z),imag(z));
        axis equal;
    end
end

function z = evaluate_interval_t(interval_list,t)
    ii_interval = 1;
    z = zeros(size(t));
    for ii = 1:length(t)
        while t(ii) > interval_list(ii_interval).t(2)
            ii_interval = ii_interval + 1;
        end
        T = interval_list(ii_interval).density_fct_handle(t(ii));
        z(ii) = interval_list(ii_interval).input_fct_handle(T);
    end
end
