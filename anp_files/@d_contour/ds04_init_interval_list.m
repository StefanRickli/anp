function [] = ds04_init_interval_list(this)
    % Inits a list of intervals that each belong to one single shape
    % function (i.e. a straight line or a circle), based on location of
    % poles and zeros on the imaginary axis and the requested radii of the
    % main half-circle and the detours.
    % 
    % PRE:  - A sorted list of purely imaginary poles and zeros.
    %       - A value for the pole and zero halfsecants that will be used.
    % 
    % POST: - A zero-initialized list of D-contour intervals.
    
    actual_halfsecant_pole =    this.halfsecant_pole;
    actual_halfsecant_zero =    this.halfsecant_zero;
    im_pz_sorted =              this.im_pz_sorted;
        
    im_poles = [im_pz_sorted([im_pz_sorted.pole] == 1).value];
    im_zeros = [im_pz_sorted([im_pz_sorted.zero] == 1).value];
    
    % This is true whenevery one of the detour circles's end points lands
    % exactly on the origin.
    pz_interval_border_on_origin = any( [(im_poles+actual_halfsecant_pole)==0,...
                                         (im_poles-actual_halfsecant_pole)==0,...
                                         (im_zeros+actual_halfsecant_zero)==0,...
                                         (im_zeros+actual_halfsecant_zero)==0]);
    
    % There are 5 base interval types along the D-contour:
    %   - a straight line from origin to the upper crop-cirlce
    %   - the upper crop-circle
    %   - the main large halfcircle
    %   - the lower crop-circle
    %   - the straight line from lower crop-circle to the origin
    % Each purely imaginary pole or zero leads to an additional straight
    % line interval (i) and its detour circle.
    % The only exeption to (i) is when a detour circle's end point lands
    % exactly on the origin.
    n_intervals = 5 + 2*length([im_poles,im_zeros]) - pz_interval_border_on_origin;
    
    this.interval_list = repmat(struct('type',[],...
                                       't',[NaN,NaN],...
                                       't_len',NaN,...
                                       'q',[NaN,NaN],...
                                       'q_len',NaN,...
                                       'density_fct_handle',[],...
                                       'input_fct_handle',[])  ,  n_intervals,1);
end