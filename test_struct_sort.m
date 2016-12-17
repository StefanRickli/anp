% Accepts two arrays containing the pole and zero values.
% Returns a struct that contains the sorted pole and zero location and
% denotes which is of which type.
function poles_zeros_sorted = test_struct_sort(poles,zeros)
    % first find all the pure imaginary poles and zeros
    poles = poles(real(poles) == 0);
    zeros = zeros(real(zeros) == 0);
    
    % strip them off the imaginary unit
    poles = imag(poles);
    zeros = imag(zeros);
    
    % remember how many purely imaginary poles and zeros we have
    npoles = length(poles);
    nzeros = length(zeros);
    
    % initialize the struct array to hold the information, courtesy to this
    % efficient way of initializing the struct array:
    % http://stackoverflow.com/questions/13664090/how-to-initialize-an-array-of-structs-in-matlab
    poles_zeros_sorted = repmat(struct('type','a','value',[]), npoles+nzeros, 1 );
    
    % fill the struct array with the pole/zero locations
    % in order to avoid for-loops, convert the num arrays first to cell
    % arrays, whose multiple outputs we then can use to fill the
    % corresponding struct array fields.
    p = repmat({'p'},npoles,1);                 % cell array full of 'p'
    pp = mat2cell(poles,1,ones(1,npoles));      % cell array containing pole values
    z = repmat({'z'},nzeros,1);                 % cell array full of 'z'
    zz = mat2cell(zeros,1,ones(1,nzeros));      % cell array containing zero values
    [poles_zeros_sorted(1:npoles).type] = p{:};
    [poles_zeros_sorted(1:npoles).value] = pp{:};
    [poles_zeros_sorted(npoles+1:end).type] = z{:};
    [poles_zeros_sorted(npoles+1:end).value] = zz{:};
    
    % Use the function below to sort the rows by the 'value'-field. Flip
    % its output in order to get a descending sorted list.
    % Now we know in which order the poles and zeros appear on the
    % imaginary axis.
    poles_zeros_sorted = flipud(sort_struct_array_by_field_index(poles_zeros_sorted,2));
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

function original_sort_example
    % http://blogs.mathworks.com/pick/2010/09/17/sorting-structure-arrays-based-on-fields/

    A = struct(...
        'name', {'mike', 'doug', 'steve', 'loren', 'jiro', 'brett', 'seth'}, ...
        'year', {2005, 2001, 1993, 1987, 2006, 2005, 1998}, ...
        'day', {'Mon', 'Fri', 'Wed', 'Fri', 'Mon', 'Mon', 'Mon'})

    tic
    Afields = fieldnames(A);
    Acell = struct2cell(A);
    sz = size(Acell)            % Notice that the this is a 3 dimensional array.
                                % For MxN structure array with P fields, the size
                                % of the converted cell array is PxMxN

    % Convert to a matrix
    Acell = reshape(Acell, sz(1), []);      % Px(MxN)

    % Make each field a column
    Acell = Acell';                         % (MxN)xP

    % Sort by first field "name"
    Acell = sortrows(Acell, 1)

    % Put back into original cell array format
    Acell = reshape(Acell', sz);

    % Convert to Struct
    Asorted = cell2struct(Acell, Afields, 1)
    toc
end