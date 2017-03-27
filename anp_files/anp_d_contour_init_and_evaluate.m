function z = anp_d_contour_init_and_evaluate(t,in_params)
    % Serves as main entry point for the creation D-contour points based on a set of input data points that are in an interval of [0,1]
    % 
    % PRE:  - t is a vector with monotonically increasing values from 0 to 1
    %       - all necessary parameters are set in in_params
    % POST: returns a vector containing the points of a close D-contour
    %       curve that begin in the origin and end in the origin
    
    % Check if 'debug_graphics' has already been set and if not, do it
    % standalone
    if any(strcmp(who('global'),'debug_graphics'))
        global debug_graphics;                                  %#ok<TLEV>
    else
        debug_graphics = false;
    end
    
    % Sort the pure imaginary poles and zeros in the order they will be
    % encountered en route of the D-contour
    in_data.im_pz_sorted = anp_in_fct_sort_pz(in_params);
    
    % Determines which combinations of pole-zero occur in which order and
    % calculates the distance between the neighbors.
    % This distances are later used to form an upper bound on the radii we
    % can use for the detours.
    in_data.im_pz_combinations = anp_in_fct_pz_combinations(in_params,in_data);
    
    [in_data.interval_list,...
     in_params.actual_separation_pole,...
     in_params.actual_separation_zero] = anp_in_fct_init_interval_list(in_params,in_data);
    
    in_data.interval_list = anp_in_fct_fill_interval_list_q_to_z(in_params,in_data);
    in_data.interval_list = anp_in_fct_fill_interval_list_t_to_q(in_params,in_data);
    
    z = evaluate_interval_t(in_data.interval_list,t);
    
    if debug_graphics
        figure;
        scatter(real(z),imag(z));
        axis equal;
    end
end

function z = evaluate_interval_t(interval_list,t)
    ii_interval = 1;
    z = zeros(size(t));
    for ii = 1:length(t)
        while (ii_interval + 1 <= length(interval_list)) && (t(ii) > interval_list(ii_interval).t(2))
            ii_interval = ii_interval + 1;
        end
        T = interval_list(ii_interval).density_fct_handle(t(ii));
        z(ii) = interval_list(ii_interval).input_fct_handle(T);
    end
end
