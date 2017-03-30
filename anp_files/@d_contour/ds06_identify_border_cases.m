function [] = ds06_identify_border_cases(this)
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
                tools.dbg('identify_border_cases: ii = %d, type = %s, value = %f\n',ii,pz_list_sorted(ii).type,pz_list_sorted(ii).value);
                error('Oops, we shouldn''t be here. Apologies! Please report this crash to ricklis@student.ethz.ch together with the input you used.');
        end
    end
end
