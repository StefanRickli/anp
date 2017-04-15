function type = get_detour_type(pz_type,remark)
    % Returns 'detour_pole' or 'detour_zero' as string or alternatively adds '_whatever-remark' to it.
    
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
            error('Oops, we shouldn''t be here. Apologies! Please report this crash to stefanrickli [at] gmx.ch together with the input you used.');
    end 
end
