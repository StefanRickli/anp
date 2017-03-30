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
