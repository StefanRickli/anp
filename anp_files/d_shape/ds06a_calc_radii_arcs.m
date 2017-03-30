function [radii,arc_lengths] = d_shape_06a_calc_radii_arcs(radii,angles,secant_pole,secant_zero)
    radii.crop = radii.inf*sin(angles.crop)/(1+sin(angles.crop));
    
    radii.detour_pole = secant_pole*sqrt(4+tan(angles.detour)^2)/4;
    radii.detour_zero = secant_zero*sqrt(4+tan(angles.detour)^2)/4;
    
    arc_lengths.inf = (pi-2*angles.crop)*radii.inf;
    arc_lengths.crop = (pi/2+angles.crop)*radii.crop;
    arc_lengths.detour_pole = 2*radii.detour_pole*asin(secant_pole/(2*radii.detour_pole));
    arc_lengths.detour_zero = 2*radii.detour_zero*asin(secant_zero/(2*radii.detour_zero));
end
