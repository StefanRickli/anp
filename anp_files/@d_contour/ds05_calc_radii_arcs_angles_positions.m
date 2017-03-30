function [] = ds05_calc_radii_arcs_angles_positions(this)
    % Calculate the parameters (radii of the detour circles from the secant lengths,
    % calculate the radii of the crop circles and their positions
    
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
