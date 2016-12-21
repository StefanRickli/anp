% muss herausfinden, wo die linearen intervalle aufhören
function interval_list = test_fill_interval_list
    clc;
    
    infinity_plus_crop_interval = 1/5;
    
    poles = 1i*[0.1,2,2.3,3,4];%,-0.24i,-0.1i,0,0.1i,0.24i,0.25i,0.255i,2i];
    zeros = 1i*[3.5];%[-2i,-0.055i,-0.05i,-0.04i,-0.01i,0,0.01i,0.04i,0.05i,0.055i,2i];
    
    pure_imag_poles = imag(poles(real(poles) == 0));
    pure_imag_zeros = imag(zeros(real(zeros) == 0));
    
    separation_max_pole = 1/4;
    separation_max_zero = 1/20;
    
    im_pz_sorted = test_struct_sort(poles,zeros);
    [interval_list,pole_zero_combinations,separation_pole,separation_zero] = test_init_interval_list(poles,zeros,im_pz_sorted,separation_max_pole,separation_max_zero);
    
%     im_pos_pz = im_pz_sorted([im_pz_sorted.value] >= 0);
%     im_neg_pz = im_pz_sorted([im_pz_sorted.value] < 0);
    
    radii.inf = 10;
    angles.crop = 7*pi/180;
    angles.detour = 45*pi/180;
    secant_pole = 2*separation_pole;
    secant_zero = 2*separation_zero;
    [radii,arc_lengths] = calculate_detour_params(radii,angles,secant_pole,secant_zero);
    positions.crop_y0 = sqrt(radii.inf^2-2*radii.crop*radii.inf);
    positions.crop_x0 = radii.crop;
    angles.detour_pole_phi0 = asin(secant_pole/(2*radii.detour_pole));
    angles.detour_zero_phi0 = asin(secant_zero/(2*radii.detour_zero));
    
    % do the first element
    % can be
    % either: a negative imaginary pz with overlap up to a positive one that
    %         just touches the origin
    % or: the first straight interval
    if ~isempty(im_pz_sorted)
        im_pz_sorted = identify_border_cases(im_pz_sorted,separation_pole,separation_zero);
    end
%     if ~isempty(im_pos_pz)
%         im_pos_pz = identify_border_cases(im_pos_pz,separation_pole,separation_zero);
%     end
%     if ~isempty(im_neg_pz)
%         im_neg_pz = identify_border_cases(im_neg_pz,separation_pole,separation_zero);
%     end
    
%     'start'
%     'end'
%     'density_fct_handle'
%     'density_fct_arguments'
%     'input_fct_handle'
%     'input_fct_arguments'
    
    interval_ii = 1;
    current_pz = NaN
    skip_one_pos_pz = false;
    pz = [im_pz_sorted.value];
    positive_idx = find(pz >= 0,1,'first');
    positive_pz_remain = length(pz)-positive_idx+1;
    prev_upper_bound = 0;
    
    if any([im_pz_sorted.neg_overlapping])
            interval_list(1).q(1) = 0;
            
            current_pz = find(pz < 0,1,'last');
            arc_length_overlapping_pole = radii.detour_pole * (asin(-im_pz_sorted(current_pz).value/radii.detour_pole) + angles.detour_pole_phi0);
            arc_length_overlapping_zero = radii.detour_zero * (asin(-im_pz_sorted(current_pz).value/radii.detour_zero) + angles.detour_zero_phi0);
            interval_length = im_pz_sorted(current_pz).pole*(arc_lengths.detour_pole - arc_length_overlapping_pole) + im_pz_sorted(current_pz).zero*(arc_lengths.detour_pole - arc_length_overlapping_zero);
            interval_list(1).q_len = interval_length;
            
            interval_list(1).q(2) = interval_length;
            prev_upper_bound = interval_list(1).q(2);
            
            switch im_pz_sorted(end).type
                case 'p'
                    interval_list(1).input_fct_handle = @(q) circ_detour(map(q,0,interval_list(1).q(2),arc_lengths.detour_pole-interval_list(1).q_len,arc_lengths.detour_pole),radii.detour_pole,secant_pole,im_pz_sorted(current_pz).value);
                case 'z'
                    interval_list(1).input_fct_handle = @(q) circ_detour(map(q,0,interval_list(1).q(2),arc_lengths.detour_zero-interval_list(1).q_len,arc_lengths.detour_zero),radii.detour_pole,secant_pole,im_pz_sorted(current_pz).value);
            end
            
            interval_ii = interval_ii + 1;
            
            
        fprintf('negative and overlapping p/z\n');
        fprintf('interval\t[%.3f\t%.3f],\tlength = %.3f,\tdetour\n',interval_list(1).q(1),interval_list(1).q(2),interval_list(1).q(2)-interval_list(1).q(1));
        
        
    elseif any([im_pz_sorted.pos_overlapping])
            current_pz = positive_idx;
            
            interval_list(1).q(1) = 0;
            
            arc_length_overlapping_pole = radii.detour_pole * (asin(-im_pz_sorted(current_pz).value/radii.detour_pole) + angles.detour_pole_phi0);
            arc_length_overlapping_zero = radii.detour_zero * (asin(-im_pz_sorted(current_pz).value/radii.detour_zero) + angles.detour_zero_phi0);
            interval_length = im_pz_sorted(current_pz).pole*(arc_lengths.detour_pole - arc_length_overlapping_pole) + im_pz_sorted(current_pz).zero*(arc_lengths.detour_pole - arc_length_overlapping_zero);
            interval_list(1).q_len = interval_length;
            
            interval_list(1).q(2) = interval_length;
            prev_upper_bound = interval_list(1).q(2);
            
            switch im_pz_sorted(end).type
                case 'p'
                    interval_list(1).input_fct_handle = @(q) circ_detour(map(q,0,interval_list(1).q(2),arc_lengths.detour_pole-interval_list(1).q_len,arc_lengths.detour_pole),radii.detour_pole,secant_pole,im_pz_sorted(current_pz).value);
                case 'z'
                    interval_list(1).input_fct_handle = @(q) circ_detour(map(q,0,interval_list(1).q(2),arc_lengths.detour_zero-interval_list(1).q_len,arc_lengths.detour_zero),radii.detour_pole,secant_pole,im_pz_sorted(current_pz).value);
            end

            positive_pz_remain = positive_pz_remain - 1;
            
            interval_ii = interval_ii + 1;
            
            
        fprintf('positive and overlapping p/z\n');
        fprintf('interval\t[%.3f\t%.3f],\tlength = %.3f,\tdetour\n',interval_list(1).q(1),interval_list(1).q(2),interval_length);
        
        
    elseif any([im_pz_sorted.pos_on_origin])
            current_pz = positive_idx;
            
            interval_list(1).q(1) = 0;
            
            interval_length = im_pz_sorted(current_pz).pole*arc_lengths.detour_pole + im_pz_sorted(current_pz).zero*arc_lengths.detour_zero;
            interval_list(1).q_len = interval_length;
            
            interval_list(1).q(2) = interval_length;
            prev_upper_bound = interval_list(1).q(2);
            
            positive_pz_remain = positive_pz_remain - 1;
            
            interval_ii = interval_ii + 1;
            
            
        fprintf('positive p/z that ends right on origin\n');
        fprintf('interval\t[%.3f\t%.3f],\tlength = %.3f,\tlinear\n',interval_list(1).q(1),interval_list(1).q(2),interval_length);
        
        
    elseif  ~isempty(positive_idx)
            prev_upper_bound = 0;
%             interval_list(1).q(1) = 0;
%             
%             interval_list(1).q_len = im_pz_sorted(positive_idx).value - separation_pole;
%             
%             interval_list(1).q(2) = im_pz_sorted(positive_idx).value - separation_pole;
%             prev_upper_bound = interval_list(1).q(2);
        fprintf('no p/z in range of origin\n');
        
        
    else
            prev_upper_bound = 0;
        
        fprintf('no imag positive p/z\n');
        
        
    end
    
    fprintf('------------------------------------------------------------------\n');
    
    % ---------------------------------------------------------------------------------------------------------------------------------------------
    % treat positive p/z

    while(positive_pz_remain)
        current_pz = current_pz + 1;
        % there are positive poles/zeros left to treat
        % make a linear region to the next one
        interval_list(interval_ii).q(1) = prev_upper_bound;
        
        if isnan(current_pz)
            current_pz = positive_idx;
            interval_length = separation_pole*im_pz_sorted(current_pz).pole + separation_zero*im_pz_sorted(current_pz).zero;
            za = 0;
        else
            interval_length = pole_zero_combinations(current_pz-1).distance - separation_pole*(im_pz_sorted(current_pz-1).pole + im_pz_sorted(current_pz).pole) - separation_zero*(im_pz_sorted(current_pz-1).zero + im_pz_sorted(current_pz).zero);            
            za = im_pz_sorted(current_pz-1).value + im_pz_sorted(current_pz-1).pole*separation_pole + im_pz_sorted(current_pz-1).zero*separation_zero;
        end
        interval_list(interval_ii).q_len = interval_length;
        
        interval_list(interval_ii).q(2) = interval_list(interval_ii).q(1) + interval_length;
        prev_upper_bound = interval_list(1).q(2);
        
        
        zb = im_pz_sorted(current_pz).value - im_pz_sorted(current_pz).pole*separation_pole + im_pz_sorted(current_pz).zero*separation_zero;
        interval_list(interval_ii).input_fct_handle = @(q) im_axis_line(q,interval_list(interval_ii).q(1),interval_list(interval_ii).q(2),za,zb);
        
        fprintf('interval\t[%.3f\t%.3f],\tlength = %.3f,\tlinear\n',interval_list(interval_ii).q(1),interval_list(interval_ii).q(2),interval_length);
        interval_ii = interval_ii + 1;
        
        % treat the pz
        interval_list(interval_ii).q(1) = prev_upper_bound;
        
        interval_length = im_pz_sorted(current_pz).pole*arc_lengths.detour_pole + im_pz_sorted(current_pz).zero*arc_lengths.detour_zero;
        interval_list(interval_ii).q_len = interval_length;
        
        interval_list(interval_ii).q(2) = interval_list(interval_ii).q(1) + interval_length;
        prev_upper_bound = interval_list(1).q(2);
        
        switch im_pz_sorted(current_pz).type
            case 'p'
                interval_list(interval_ii).input_fct_handle = @(q) circ_detour(map(q,interval_list(interval_ii).q(1),interval_list(interval_ii).q(2),0,arc_lengths.detour_pole),radii.detour_pole,secant_pole,im_pz_sorted(current_pz).value);
            case 'z'
                interval_list(interval_ii).input_fct_handle = @(q) circ_detour(map(q,interval_list(interval_ii).q(1),interval_list(interval_ii).q(2),0,arc_lengths.detour_zero),radii.detour_pole,secant_pole,im_pz_sorted(current_pz).value);
        end
        
        fprintf('interval\t[%.3f\t%.3f],\tlength = %.3f,\tdetour\n',interval_list(interval_ii).q(1),interval_list(interval_ii).q(2),interval_length);
        interval_ii = interval_ii + 1;
        
        positive_pz_remain = positive_pz_remain - 1;
    end
    
    % do the last linear interval before crop1
    interval_list(interval_ii).q(1) = prev_upper_bound;

    interval_length = positions.crop_y0 - (im_pz_sorted(current_pz).value + im_pz_sorted(current_pz).pole*separation_pole + im_pz_sorted(current_pz).zero*separation_zero);
    interval_list(interval_ii).q_len = interval_length;

    interval_list(interval_ii).q(2) = interval_list(interval_ii).q(1) + interval_length;

    za = im_pz_sorted(current_pz).value + im_pz_sorted(current_pz).pole*separation_pole + im_pz_sorted(current_pz).zero*separation_zero;
    zb = positions.crop_y0;
    interval_list(interval_ii).input_fct_handle = @(q) im_axis_line(q,interval_list(interval_ii).q(1),interval_list(interval_ii).q(2),za,zb);

    fprintf('interval\t[%.3f\t%.3f],\tlength = %.3f,\tlinear_last_pos\n',interval_list(interval_ii).q(1),interval_list(interval_ii).q(2),interval_length);
    interval_ii = interval_ii + 1;
    
    
    % crop1
    interval_list(interval_ii).q(1) = interval_list(interval_ii-1).q(2);
    interval_list(interval_ii).q_len = arc_lengths.crop;
    interval_list(interval_ii).q(2) = interval_list(interval_ii).q(1) + arc_lengths.crop;
    
    interval_list(interval_ii).input_fct_handle = @(q) circ_normal(-map(q,interval_list(interval_ii).q(1),interval_list(interval_ii).q(2),0,arc_lengths.crop),radii.crop,pi,positions.crop_x0,positions.crop_y0);
    
    fprintf('interval\t[%.3f\t%.3f],\tlength = %.3f,\tcrop\n',interval_list(interval_ii).q(1),interval_list(interval_ii).q(2),arc_lengths.crop);
    interval_ii = interval_ii + 1;
    
    % inf
    interval_list(interval_ii).q(1) = interval_list(interval_ii-1).q(2);
    interval_list(interval_ii).q_len = arc_lengths.inf;
    interval_list(interval_ii).q(2) = interval_list(interval_ii).q(1) + arc_lengths.inf;
    
    interval_list(interval_ii).input_fct_handle = @(q) circ_normal(-map(q,interval_list(interval_ii).q(1),interval_list(interval_ii).q(2),0,arc_lengths.inf),radii.inf,pi/2-angles.crop,0,0);    
    
    fprintf('interval\t[%.3f\t%.3f],\tlength = %.3f,\tinf\n',interval_list(interval_ii).q(1),interval_list(interval_ii).q(2),arc_lengths.inf);
    interval_ii = interval_ii + 1;
    
    % crop2
    interval_list(interval_ii).q(1) = interval_list(interval_ii-1).q(2);
    interval_list(interval_ii).q_len = arc_lengths.crop;
    interval_list(interval_ii).q(2) = interval_list(interval_ii).q(1) + arc_lengths.crop;
    
    interval_list(interval_ii).input_fct_handle = @(q) circ_normal(-map(q,interval_list(interval_ii).q(1),interval_list(interval_ii).q(2),0,arc_lengths.crop),radii.crop,-pi/2+angles.crop,positions.crop_x0,-positions.crop_y0);
    
    fprintf('interval\t[%.3f\t%.3f],\tlength = %.3f,\tcrop\n',interval_list(interval_ii).q(1),interval_list(interval_ii).q(2),arc_lengths.crop);
    interval_ii = interval_ii + 1;
    
    % treat negative p/z

    % then: for-loop through every element with lower index than the border
    % case
    
    if any([im_pz_sorted.pos_overlapping])
        fprintf('treating a positive and overlapping p/z\n');
        
    elseif any([im_pz_sorted.neg_on_origin])
        fprintf('treating a negative p/z that ends right on origin\n');
        
    end
    
    
    for kk = 1:15
        t = intvl(interval_list(kk).q);
        plot(interval_list(kk).input_fct_handle(t));
        hold on;
    end
    
    axis equal;
    
end

function pz_list_sorted = identify_border_cases(pz_list_sorted,separation_pole,separation_zero)
    for ii = 1:length(pz_list_sorted)
        pz_value = pz_list_sorted(ii).value;
        switch(pz_list_sorted(ii).type)
            case 'p'
                if pz_value < 0 && (pz_value + separation_pole == 0)
                    pz_list_sorted(ii).neg_on_origin = true;
                elseif pz_value < 0 && (pz_value + separation_pole > 0)
                    pz_list_sorted(ii).neg_overlapping = true;
                elseif pz_value >= 0 && (pz_value - separation_pole == 0)
                    pz_list_sorted(ii).pos_on_origin = true;
                elseif pz_value >= 0 && (pz_value - separation_pole < 0)
                    pz_list_sorted(ii).pos_overlapping = true;
                end
            case 'z'
                if pz_value < 0 && (pz_value + separation_zero == 0)
                    pz_list_sorted(ii).neg_on_origin = true;
                elseif pz_value < 0 && (pz_value + separation_zero > 0)
                    pz_list_sorted(ii).neg_overlapping = true;
                elseif pz_value >= 0 && (pz_value - separation_zero == 0)
                    pz_list_sorted(ii).pos_on_origin = true;
                elseif pz_value >= 0 && (pz_value - separation_zero < 0)
                    pz_list_sorted(ii).pos_overlapping = true;
                end
            otherwise
                fprintf('identify_border_cases: ii = %d, type = %s, value = %f\n',ii,pz_list_sorted(ii).type,pz_list_sorted(ii).value);
                error('Oops, we shouldn''t be here. Apologies! Please report this crash to ricklis@student.ethz.ch together with the input you used.');
        end
    end
end

function [radii,arc_lengths] = calculate_detour_params(radii,angles,secant_pole,secant_zero)
    radii.crop = radii.inf*sin(angles.crop)/(1+sin(angles.crop));
    
    radii.detour_pole = secant_pole*sqrt(4+tan(angles.detour)^2)/4;
    radii.detour_zero = secant_zero*sqrt(4+tan(angles.detour)^2)/4;
    
    arc_lengths.inf = (pi-2*angles.crop)*radii.inf;
    arc_lengths.crop = (pi/2+angles.crop)*radii.crop;
    arc_lengths.detour_pole = 2*radii.detour_pole*asin(secant_pole/(2*radii.detour_pole));
    arc_lengths.detour_zero = 2*radii.detour_zero*asin(secant_zero/(2*radii.detour_zero));
end

function z = circ_normal(q,R,phi_0,x0,y0)
    assert((real(y0) == 0 && imag(y0) ~= 0) || (real(y0) ~= 0 && imag(y0) == 0) || (y0 == 0));
    y0 = real(y0) + imag(y0);
    
    z = R*exp(1i*(phi_0 + q/R)) + x0 + 1i*y0;
end

function z = circ_detour(q,R,sec,y_pz)
    assert((real(y_pz) == 0 && imag(y_pz) ~= 0) || (real(y_pz) ~= 0 && imag(y_pz) == 0));
    y_pz = real(y_pz) + imag(y_pz);
    
    z = R*(exp(1i*(q/R - asin(sec/(2*R)))) - sqrt(1-(sec/(2*R))^2)) + 1i*y_pz;
end

function z = im_axis_line(q,qa,qb,za,zb)
    z = 1i*map(q,qa,qb,za,zb);
end

function y = map(x,t0,t1,u0,u1)
    y = ((u0-u1).*x + (t0*u1-t1*u0))/(t0-t1);
end

function t = intvl(in)
    t = linspace(in(1),in(2),100);
end
