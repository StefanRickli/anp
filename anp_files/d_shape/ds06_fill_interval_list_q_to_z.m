function interval_list = d_shape_06_fill_interval_list_q_to_z(in_params,in_data)
    
    if any(strcmp(who('global'),'debug_graphics'))
        global debug_graphics;
    else
        debug_graphics = false;
    end
    
    % Copy into local variables to shorten the code
    radii =                     in_params.radii;
    angles =                    in_params.angles;
    halfsecant_pole =           in_params.actual_halfsecant_pole;
    halfsecant_zero =           in_params.actual_halfsecant_zero;    
    im_pz_sorted =              in_data.im_pz_sorted;
    pole_zero_combinations =    in_data.im_pz_combinations;
    interval_list =             in_data.interval_list;
    
    % Calculate the parameters (radii of the detour circles from the secant lengths,
    % calculate the radii of the crop circles and their positions
    secant_pole =               2 * halfsecant_pole;
    secant_zero =               2 * halfsecant_zero;
    [radii,arc_lengths] =       d_shape_06a_calc_radii_arcs(radii,angles,secant_pole,secant_zero);
    angles.detour_pole_phi0 =   asin(secant_pole / (2*radii.detour_pole));
    angles.detour_zero_phi0 =   asin(secant_zero / (2*radii.detour_zero));
    positions.crop_y0 =         sqrt(radii.inf^2 - 2*radii.crop*radii.inf);
    positions.crop_x0 =         radii.crop;
    
    if ~isempty(im_pz_sorted)
        im_pz_sorted = identify_border_cases(im_pz_sorted,halfsecant_pole,halfsecant_zero);
    end
    
    d_shape_06a_1st_interval();
    d_shape_06b_positive_pz();
    d_shape_06c_last_straight_before_crop1();
    d_shape_06d_crop1_inf_crop2();
    d_shape_06g_first_after_crop2();
    d_shape_06h_negative_pz();
    d_shape_06i_last_interval();
    
    if debug_graphics
        figure;
        for kk = 1:length(interval_list)
            t = intvl(interval_list(kk).q);
            scatter(real(interval_list(kk).input_fct_handle(t)),imag(interval_list(kk).input_fct_handle(t)));
            hold on;
        end
        axis equal;
    end
end

function pz_list_sorted = identify_border_cases(pz_list_sorted,halfsecant_pole,halfsecant_zero)
    for ii = 1:length(pz_list_sorted)
        pz_value = pz_list_sorted(ii).value;
        switch(pz_list_sorted(ii).type)
            case 'p'
                if pz_value < 0 && (pz_value + halfsecant_pole == 0)
                    pz_list_sorted(ii).neg_on_origin = true;
                elseif pz_value < 0 && (pz_value + halfsecant_pole > 0)
                    pz_list_sorted(ii).neg_overlapping = true;
                elseif pz_value >= 0 && (pz_value - halfsecant_pole == 0)
                    pz_list_sorted(ii).pos_on_origin = true;
                elseif pz_value >= 0 && (pz_value - halfsecant_pole < 0)
                    pz_list_sorted(ii).pos_overlapping = true;
                end
            case 'z'
                if pz_value < 0 && (pz_value + halfsecant_zero == 0)
                    pz_list_sorted(ii).neg_on_origin = true;
                elseif pz_value < 0 && (pz_value + halfsecant_zero > 0)
                    pz_list_sorted(ii).neg_overlapping = true;
                elseif pz_value >= 0 && (pz_value - halfsecant_zero == 0)
                    pz_list_sorted(ii).pos_on_origin = true;
                elseif pz_value >= 0 && (pz_value - halfsecant_zero < 0)
                    pz_list_sorted(ii).pos_overlapping = true;
                end
            otherwise
                tools.dbg('identify_border_cases: ii = %d, type = %s, value = %f\n',ii,pz_list_sorted(ii).type,pz_list_sorted(ii).value);
                error('Oops, we shouldn''t be here. Apologies! Please report this crash to ricklis@student.ethz.ch together with the input you used.');
        end
    end
end

function type = get_detour_type(pz_type,remark)
    underline = [];
    if ~isempty(remark)
        underline = '_';
    end
    
    switch pz_type
        case 'p'
            type = ['detour_pole',underline,remark];
        case 'z'
            type = ['detour_zero',underline,remark];
        otherwise
            error('Oops, we shouldn''t be here. Apologies! Please report this crash to ricklis@student.ethz.ch together with the input you used.');
    end 
end

function z = circ_normal(q,R,phi_0,x0,y0)
    assert((real(y0) == 0 && imag(y0) ~= 0) || (real(y0) ~= 0 && imag(y0) == 0) || (y0 == 0));
    y0 = real(y0) + imag(y0);
    
    z = R*exp(1i*(phi_0 + q/R)) + x0 + 1i*y0;
end

function z = circ_detour(q,R,sec,y_pz)
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