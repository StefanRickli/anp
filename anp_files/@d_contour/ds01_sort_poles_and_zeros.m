% Accepts two arrays containing the pole and zero values.
% Returns a struct that contains the sorted pole and zero location and
% denotes which is of which type.
function [] = ds01_sort_poles_and_zeros(this)
    poles = this.poles;
    zeros = this.zeros;
    
    % first find all the pure imaginary poles and zeros
    tol =   100*eps;
    poles = poles(abs(real(poles)) < tol);
    zeros = zeros(abs(real(zeros)) < tol);
    
    % strip them off the imaginary unit
    poles = imag(poles);
    zeros = imag(zeros);
    
    % remember how many purely imaginary poles and zeros we have
    npoles = length(poles);
    nzeros = length(zeros);
    
    if npoles + nzeros == 0
        this.im_pz_sorted = struct('type',{},'pole',{},'zero',{},'value',{},'neg_on_origin',{},'pos_on_origin',{},'neg_overlapping',{},'pos_overlapping',{});
        return;
    end
    
    % initialize the struct array to hold the information, courtesy to this
    % efficient way of initializing the struct array:
    % http://stackoverflow.com/questions/13664090/how-to-initialize-an-array-of-structs-in-matlab
    % 
    % Meaning of the fields:
    % type:
    %   (p)ole, (z)ero; it's there for the debugging human
    % pole/zero:
    %   booleans representing the type of point
    % neg_overlapping, pos_overlapping:
    %   These fields will be set later. They indicate the special case when
    %   an imaginary p/z is so close to the origin such that the detour
    %   circle gets split into a part that resid in the lower and upper
    %   half plane.
    % neg_on_origin, pos_on_origin:
    %   Similar to the former two fields these booleans indicate the case
    %   when an imaginary p/z detour circle's end point lies exactly on the
    %   origin. In this case it doesn't get split in two parts and there
    %   will be no straight line from the end point of the detour circle to
    %   the origin.
    poles_zeros_sorted = repmat(struct('type','a','pole',false,'zero',false,'value',[],'neg_on_origin',false,'pos_on_origin',false,'neg_overlapping',false,'pos_overlapping',false), npoles+nzeros, 1 );
    
    % First fill the struct array 'poles_zeros_sorted' with the pole/zero
    % locations, unsorted, first poles then zeros
    % 
    % (In order to avoid for-loops, convert the num arrays first to cell
    %  arrays, whose multiple outputs we then can use to fill the
    %  corresponding struct array fields.)
    if npoles > 0
        type =      repmat({'p'},1,npoles);                 % 1D cell array full of 'p'
        pole =      repmat({true},1,npoles);                % 1D cell array full of boolean true
        value =     mat2cell(poles,1,ones(1,npoles));       % 1D cell array containing the pole values
        
        % Fill the struct array with the previously prepared data
        [poles_zeros_sorted(1:npoles).type] = type{:};
        [poles_zeros_sorted(1:npoles).pole] = pole{:};
        [poles_zeros_sorted(1:npoles).value] = value{:};
    end
    if nzeros > 0
        type =      repmat({'z'},1,nzeros);                 % cell array full of 'z'
        zero =      repmat({true},1,nzeros);                 % cell array full of 'z'
        value =     mat2cell(zeros,1,ones(1,nzeros));      % cell array containing zero values
        
        % Amend the struct array with the previously prepared data
        [poles_zeros_sorted(npoles+1:end).type] = type{:};
        [poles_zeros_sorted(npoles+1:end).zero] = zero{:};
        [poles_zeros_sorted(npoles+1:end).value] = value{:};
    end
    
    % Use the function below to sort the rows by the 'value'-field.
    % Now we know in which order the poles and zeros appear on the
    % imaginary axis.
    if npoles+nzeros > 1
        poles_zeros_sorted = sort_struct_array_by_field_index(poles_zeros_sorted,4);
        poles_zeros_sorted = check_for_multiple_pz_at_same_value(poles_zeros_sorted);
    end
    
    this.im_pz_sorted = poles_zeros_sorted;
end

function sorted_struct_array = sort_struct_array_by_field_index(array_in,ind)
    % Code adopted from
    % http://blogs.mathworks.com/pick/2010/09/17/sorting-structure-arrays-based-on-fields/
    
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
    % Removes any poles or zeros that occur more than once at the exact
    % same spot/value.
    % Note: We don't identify for pole/zero cancellations. If one happens,
    %       the one to appear later in the list will get dropped and we
    %       still make a detour.
    
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