function [] = ds05_calc_radii_arcs_angles_positions(this)
    % Calculates the parameters (radii of the detour circles from the secant lengths, radii of the crop circles and their positions).
    
    %  ^ Im{}
    %  |
    %  |
    %  |         ----\
    %  |      -     /   ---\
    %  |  /  crop  /         --\
    %  /---------c               --\
    %  |                             -\
    %  |                                \
    %  |                                 \
    %  |-                                 |
    %  |  \                                |
    % p+---|  pole detour                   |
    %  |  /                                 |
    %  |-                                   |
    %  |                                     |
    %  |-                                    |
    %  |  \                                   |
    % z+---| zero detour                      |
    %  |  /                                   |  main (~inf) halfcircle
    %  |-                                     |
    %  |                                      |
    %  m-------------------------------------------------------------> Re{}
    % 
    % Detour arcs:
    % Be aware that the center of the detour arcs are NOT on the Im-axis,
    % but somewhat on the LHP side of it. This allows for a non-90° (and
    % thus smoother) change of direction when we change from a chunk of
    % straight axis to the detour.
    % This however leads to a slightly more complicated calculation of the
    % parametrization of those detours via their secants that lie exactly
    % on the Im-axis.
    % 
    % Crop arcs:
    % The purpose of the crop arcs is to smoothen the change of direction
    % from the last straight axis part to the main halfcircle. Their
    % parameters are chosen such that the main halfcircle loses angles.crop
    % degrees to the crop arcs at both of its ends.
    
    this.secant_pole =             	2 * this.halfsecant_pole;
    this.secant_zero =          	2 * this.halfsecant_zero;
    
    this.radii.crop =               this.radii.inf * sin(this.angles.crop) / (1 + sin(this.angles.crop));
    
    this.radii.detour_pole =        this.secant_pole * sqrt(4 + tan(this.angles.detour)^2 ) / 4;
    this.radii.detour_zero =        this.secant_zero * sqrt(4 + tan(this.angles.detour)^2 ) / 4;
    
    this.arc_lengths.inf =          (pi - 2*this.angles.crop) * this.radii.inf;
    this.arc_lengths.crop =         (pi/2 + this.angles.crop) * this.radii.crop;
    this.arc_lengths.detour_pole =  2 * this.radii.detour_pole * asin(this.secant_pole / (2 * this.radii.detour_pole) );
    this.arc_lengths.detour_zero =  2 * this.radii.detour_zero * asin(this.secant_zero / (2 * this.radii.detour_zero) );
    
    this.angles.detour_pole_phi0 =  asin(this.secant_pole / (2 * this.radii.detour_pole) );
    this.angles.detour_zero_phi0 =  asin(this.secant_zero / (2 * this.radii.detour_zero) );
    
    this.positions.crop_y0 =        sqrt( (this.radii.inf)^2 - (2 * this.radii.crop * this.radii.inf) );
    this.positions.crop_x0 =        this.radii.crop;
end
