
poles = [10,3,-1];
zeros = [2,-5];

npoles = length(poles);
nzeros = length(zeros);

% http://stackoverflow.com/questions/13664090/how-to-initialize-an-array-of-structs-in-matlab
poles_zeros_sorted = repmat(struct('type','a','value',[]), npoles+nzeros, 1 );

p = repmat({'p'},npoles,1);
pp = mat2cell(poles,1,ones(1,npoles));
z = repmat({'z'},nzeros,1);
zz = mat2cell(zeros,1,ones(1,nzeros));
[poles_zeros_sorted(1:npoles).type] = p{:};
[poles_zeros_sorted(1:npoles).value] = pp{:};
[poles_zeros_sorted(npoles+1:end).type] = z{:};
[poles_zeros_sorted(npoles+1:end).value] = zz{:};

poles_zeros_sorted = flipud(sort_struct_array_by_fieldind(poles_zeros_sorted,2));

% code adopted from
% http://blogs.mathworks.com/pick/2010/09/17/sorting-structure-arrays-based-on-fields/
function sorted_struct_array = sort_struct_array_by_fieldind(array_in,ind)
    Afields = fieldnames(array_in);
    
    if ind > length(Afields)
        error('You used an index larger than the number of fields of this struct, noob! :-)');
    end
    
    Acell = struct2cell(array_in);
    sz = size(Acell)            % Notice that the this is a 3 dimensional array.
                                % For MxN structure array with P fields, the size
                                % of the converted cell array is PxMxN

    % Convert to a matrix
    Acell = reshape(Acell, sz(1), []);      % Px(MxN)

    % Make each field a column
    Acell = Acell';                         % (MxN)xP

    % Sort by field #ind
    Acell = sortrows(Acell, ind)

    % Put back into original cell array format
    Acell = reshape(Acell', sz);

    % Convert to Struct
    sorted_struct_array = cell2struct(Acell, Afields, 1)

end

function original_example
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