function [interval_ii] = ds07d_crop1_inf_crop2(this,interval_ii)
    % Handles the two crop arcs and the main 'inf' halfcircle.
    %  ^ Im{}
    %  |
    %  |
    %  |         ----\
    %  |      -     /   ---\
    %  |  /  crop1 /         --\
    %  /----------               --\
    %  |                             -\
    %  |                                \
    %  |                                 \
    %  |                                  |
    %  |                                   |
    %  |                                    |
    %  |                                     |
    %  |                                     |
    %  |                                     |  main (~inf) ~halfcircle
    % -+-------------------------------------------------------------> Re{}
    %  |                                     |
    %  |                                     |
    %  |                                     |
    %  |                                    |
    %  |                                   |
    %  |                                  |
    %  |                                 /
    %  |                                /
    %  |                             -/
    %  \----------               --/
    %  |  \  crop2 \         --/
    %  |      -     \   ---/
    %  |         ----/
    %  |
    
    %% Init
    
    % Make the used variables local to get rid of 'this.'
    interval_list =         this.interval_list;
    arc_lengths =           this.arc_lengths;
    radii =                 this.radii;
    positions =             this.positions;
    angles =                this.angles;
    
    
    %% Crop 1
    
    
    interval_list(interval_ii).type = 'crop';
    
    %           q(1)   q_len   q(2)
    % ------------|--------------|-------------> q
    % prev_intvl    current_intvl  next_intvl
    interval_list(interval_ii).q(1) = interval_list(interval_ii-1).q(2);
    interval_list(interval_ii).q_len = arc_lengths.crop;
    interval_list(interval_ii).q(2) = interval_list(interval_ii).q(1) + arc_lengths.crop;
    
    % Construct a function handle which will be fed with q = 0 to
    % q = 'length of this interval'.
    % Here we predefine a regular complex exponential function with the
    % necessary arguments like radius, initial angle and center position.
    interval_list(interval_ii).input_fct_handle = @(q) circ_normal(-map(q,interval_list(interval_ii).q(1),interval_list(interval_ii).q(2),0,arc_lengths.crop), ...
                                                                   radii.crop, ...
                                                                   pi, ...
                                                                   positions.crop_x0, ...
                                                                   positions.crop_y0);
    
    dbg_out('interval\t[%.3f\t%.3f],\tlength = %.3f,\tcrop1\n',interval_list(interval_ii).q(1),interval_list(interval_ii).q(2),arc_lengths.crop);
    interval_ii = interval_ii + 1;
    
    
    %% Main (inf) halfcircle
    
    % Refer to documentation of Crop 1 interval
    interval_list(interval_ii).type = 'inf';
    interval_list(interval_ii).q(1) = interval_list(interval_ii-1).q(2);
    interval_list(interval_ii).q_len = arc_lengths.inf;
    interval_list(interval_ii).q(2) = interval_list(interval_ii).q(1) + arc_lengths.inf;
    
    interval_list(interval_ii).input_fct_handle = @(q) circ_normal(-map(q,interval_list(interval_ii).q(1),interval_list(interval_ii).q(2),0,arc_lengths.inf), ...
                                                                   radii.inf, ...
                                                                   pi/2-angles.crop, ...
                                                                   0, ...
                                                                   0);
    
    dbg_out('interval\t[%.3f\t%.3f],\tlength = %.3f,\tinf\n',interval_list(interval_ii).q(1),interval_list(interval_ii).q(2),arc_lengths.inf);
    interval_ii = interval_ii + 1;
    
    
    %% Crop 2
    
    % Refer to documentation of Crop 1 interval
    interval_list(interval_ii).type = 'crop';
    interval_list(interval_ii).q(1) = interval_list(interval_ii-1).q(2);
    interval_list(interval_ii).q_len = arc_lengths.crop;
    interval_list(interval_ii).q(2) = interval_list(interval_ii).q(1) + arc_lengths.crop;
    
    interval_list(interval_ii).input_fct_handle = @(q) circ_normal(-map(q,interval_list(interval_ii).q(1),interval_list(interval_ii).q(2),0,arc_lengths.crop), ...
                                                                   radii.crop, ...
                                                                   angles.crop - pi/2, ...
                                                                    positions.crop_x0, ...
                                                                   -positions.crop_y0);
    
    dbg_out('interval\t[%.3f\t%.3f],\tlength = %.3f,\tcrop2\n',interval_list(interval_ii).q(1),interval_list(interval_ii).q(2),arc_lengths.crop);
    interval_ii = interval_ii + 1;
    
    
    %% Cleanup
    
    % Write back changes to the interval list.
    this.interval_list = interval_list;
end