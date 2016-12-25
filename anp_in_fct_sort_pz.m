% Accepts two arrays containing the pole and zero values.
% Returns a struct that contains the sorted pole and zero location and
% denotes which is of which type.
function poles_zeros_sorted = anp_in_fct_sort_pz(in_params)
    poles = in_params.poles;
    zeros = in_params.zeros;
    
    % first find all the pure imaginary poles and zeros
    poles = poles(real(poles) == 0);
    zeros = zeros(real(zeros) == 0);
    
    % strip them off the imaginary unit
    poles = imag(poles);
    zeros = imag(zeros);
    
    % remember how many purely imaginary poles and zeros we have
    npoles = length(poles);
    nzeros = length(zeros);
    
    if npoles + nzeros == 0
        poles_zeros_sorted = struct('type',{},'pole',{},'zero',{},'value',{},'neg_on_origin',{},'pos_on_origin',{},'neg_overlapping',{},'pos_overlapping',{});
        return;
    end
    
    % initialize the struct array to hold the information, courtesy to this
    % efficient way of initializing the struct array:
    % http://stackoverflow.com/questions/13664090/how-to-initialize-an-array-of-structs-in-matlab
    poles_zeros_sorted = repmat(struct('type','a','pole',false,'zero',false,'value',[],'neg_on_origin',false,'pos_on_origin',false,'neg_overlapping',false,'pos_overlapping',false), npoles+nzeros, 1 );
    
    % fill the struct array with the pole/zero locations
    % in order to avoid for-loops, convert the num arrays first to cell
    % arrays, whose multiple outputs we then can use to fill the
    % corresponding struct array fields.
    if npoles > 0
        p = repmat({'p'},1,npoles);                 % cell array full of 'p'
        pp = repmat({true},1,npoles);
        ppp = mat2cell(poles,1,ones(1,npoles));      % cell array containing pole values
        [poles_zeros_sorted(1:npoles).type] = p{:};
        [poles_zeros_sorted(1:npoles).pole] = pp{:};
        [poles_zeros_sorted(1:npoles).value] = ppp{:};
    end
    if nzeros > 0
        z = repmat({'z'},1,nzeros);                 % cell array full of 'z'
        zz = repmat({true},1,nzeros);                 % cell array full of 'z'
        zzz = mat2cell(zeros,1,ones(1,nzeros));      % cell array containing zero values
        [poles_zeros_sorted(npoles+1:end).type] = z{:};
        [poles_zeros_sorted(npoles+1:end).zero] = zz{:};
        [poles_zeros_sorted(npoles+1:end).value] = zzz{:};
    end
    
    % Use the function below to sort the rows by the 'value'-field.
    % Now we know in which order the poles and zeros appear on the
    % imaginary axis.
    if npoles+nzeros > 1
        poles_zeros_sorted = sort_struct_array_by_field_index(poles_zeros_sorted,4);
        poles_zeros_sorted = check_for_multiple_pz_at_same_value(poles_zeros_sorted);
    end
end

% code adopted from
% http://blogs.mathworks.com/pick/2010/09/17/sorting-structure-arrays-based-on-fields/
function sorted_struct_array = sort_struct_array_by_field_index(array_in,ind)
    Afields = fieldnames(array_in);
    
    if ind > length(Afields)
        error('You used an index larger than the number of fields of this struct, noob! :-)');
    end
    
    Acell = struct2cell(array_in);
    sz = size(Acell);           % Notice that the this is a 3 dimensional array.
                                % For MxN structure array with P fields, the size
                                % of the converted cell array is PxMxN

    % Convert to a matrix
    Acell = reshape(Acell, sz(1), []);      % Px(MxN)

    % Make each field a column
    Acell = Acell';                         % (MxN)xP

    % Sort by field #ind
    Acell = sortrows(Acell, ind);

    % Put back into original cell array format
    Acell = reshape(Acell', sz);

    % Convert to Struct
    sorted_struct_array = cell2struct(Acell, Afields, 1);

end

function poles_zeros_sorted = check_for_multiple_pz_at_same_value(poles_zeros_sorted)
    start_over = false;
    while(true)
        for ii = 2:length(poles_zeros_sorted)
            if poles_zeros_sorted(ii-1).value == poles_zeros_sorted(ii).value
                poles_zeros_sorted(ii-1) = [];
                start_over = true;
                break;
            end
        end
        if start_over
            start_over = false;
            continue;
        end
        break;
    end
end