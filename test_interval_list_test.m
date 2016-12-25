function test_interval_list_test
    clc;
    global debug
    debug = false;
    
    poles = 1i*[];
    zeros = 1i*[];
    %poles = 1i*[0,7];
    %zeros = 1i*[-5];
    radii.inf = 10;
    angles.crop = 7*pi/180;
    angles.detour = 45*pi/180;    
    separation_max_pole = 1/4;
    separation_max_zero = 1/20;

    interval_list = test_fill_interval_list_t_to_q(poles,zeros,radii,angles,separation_max_pole,separation_max_zero);
    
%     q_max = interval_list(end).q(2);
%     t = 0:1/1000:1;
%     q = map(t,0,1,0,q_max);
%     z = evaluate_interval_q(interval_list,q);
    t = 0:1/1000:1;
    z = evaluate_interval_t(interval_list,t(1:end-1));
    
    figure;
    scatter(real(z),imag(z));
    axis equal;
end

function z = evaluate_interval_t(interval_list,t)
    ii_interval = 1;
    z = zeros(size(t));
    for ii = 1:length(t)
        while t(ii) > interval_list(ii_interval).t(2)
            ii_interval = ii_interval + 1;
        end
%         ta = interval_list(ii_interval).t(1);
%         tb = interval_list(ii_interval).t(2);
%         qa = interval_list(ii_interval).q(1);
%         qb = interval_list(ii_interval).q(2);
        T = interval_list(ii_interval).density_fct_handle(t(ii));
        z(ii) = interval_list(ii_interval).input_fct_handle(T);
    end
end


function z = evaluate_interval_q(interval_list,q)
    ii_interval = 1;
    z = zeros(size(q));
    for ii = 1:length(q)
        while q(ii) > interval_list(ii_interval).q(2)
            ii_interval = ii_interval + 1;
        end
        z(ii) = interval_list(ii_interval).input_fct_handle(q(ii));
    end
end

function y = map(x,t0,t1,u0,u1)
    y = ((u0-u1).*x + (t0*u1-t1*u0))/(t0-t1);
end

