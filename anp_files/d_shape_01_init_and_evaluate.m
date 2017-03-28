function z = d_shape_01_init_and_evaluate(t,in_params)
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
    in_data.im_pz_sorted = d_shape_02_sort_pz(in_params);
    
    % Determines which combinations of pole-zero occur in which order and
    % calculates the distance between the neighbors.
    % These distances are later used to form an upper bound on the radii we
    % can use for the detours.
    in_data.im_pz_combinations = d_shape_03_pz_combinations(in_params,in_data);
    
    % Based on location of poles and zeros on the imaginary axis and the
    % requested radii of the main half-circle and the detours create a list
    % of intervals that will segment [0,1] into blocks that each belong to
    % one single shape function, i.e. a straight line or a circle.
    % 
    % Also as we enforce a piece of straight line between neighboring
    % detours, the detour secants might be upper bounded by the minimum
    % separation of imaginary poles/zeros. We're interested in finding out
    % the maximal allowed detour secant which will be defined in terms of
    % the 'actual_halfsecant_...' variables later.
    [in_data.interval_list,...
     in_params.actual_halfsecant_pole,...
     in_params.actual_halfsecant_zero] = d_shape_04_init_interval_list(in_params,in_data);
    
    % First we create appropriate functions that map an input to a certain
    % shape, i.e. a straight line or a circle beginning at theta degrees,
    % and fill them into a handle list (functions q |--> z).
    % Then we stack up the length that the plot lines of those functions
    % will take (for straight lines just take the difference between
    % starting and end point, for circles take their arc length) from zero
    % up to the total length of the D-contour shape. This list of intervals
    % will be the input to the functions q |--> z.
    in_data.interval_list = d_shape_05_fill_interval_list_q_to_z(in_params,in_data);
    
    % As tf_processor specifies inputs in an interval of t in [0,1], we
    % need to map those values to the q-values. Do this such that many
    % t-values lie inside q-intervals that cause much movement in the
    % z-values, i.e. scale the spatial resolution.
    % We do this by specifying appropriate point density functions that
    % map parts of [0,1] nonlinearly to the q-intervals
    % (functions t |--> q).
    in_data.interval_list = d_shape_06_fill_interval_list_t_to_q(in_params,in_data);
    
    % Finally iterate through every data point that we got from
    % tf_processor and calculate its corresponding value in the z-plane by
    % performing the mapping t |--> q |--> z in the appropriate interval.
    z = evaluate_interval_t(in_data.interval_list,t);
    
    if debug_graphics
        figure;
        scatter(real(z),imag(z));
        axis equal;
    end
end

function z = evaluate_interval_t(interval_list,t)
    % Performs the mapping t |--> q |--> z for every point in the vector t
    % 
    % PRE: Monotonically increasing list of real values in the range [0,1]
    
    % As we know that the values in t are monotonically increasing, there's
    % no need to look for the right interval for every point of t. Simply
    % check whether t(ii) > current interval's upper limit and switch to
    % the next interval if the comparison yields true.
    ii_interval = 1;
    z = zeros(size(t));
    for ii = 1:length(t)
        while (ii_interval + 1 <= length(interval_list)) && (t(ii) > interval_list(ii_interval).t(2))
            ii_interval = ii_interval + 1;
        end
        
        % Map t(ii) |--> q(ii)
        Q = interval_list(ii_interval).density_fct_handle(t(ii));
        
        % Map q(ii) |--> z(ii)
        z(ii) = interval_list(ii_interval).input_fct_handle(Q);
    end
end
