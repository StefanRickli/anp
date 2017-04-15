function [] = ds06_identify_border_cases(this)
    % Checks for every purely imaginary p/z in the list if its detour arc coincidentally crosses the real axis or lands exactly on it.
    % We need this info later to plan the intervals of the D-contour.
    
    if isempty(this.im_pz_sorted)
        return;
    end
    
    pz_list_sorted =    this.im_pz_sorted;
    halfsecant_pole =   this.halfsecant_pole;
    halfsecant_zero =   this.halfsecant_zero;
    
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
                dbg_out('identify_border_cases: ii = %d, type = %s, value = %f\n',ii,pz_list_sorted(ii).type,pz_list_sorted(ii).value);
                error('Oops, we shouldn''t be here. Apologies! Please report this crash to stefanrickli [at] gmx.ch together with the input you used.');
        end
    end
    
    this.im_pz_sorted = pz_list_sorted;
end
