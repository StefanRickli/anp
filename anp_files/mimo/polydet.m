
% The contents of this cell array represent polynomial coefficients.
A = {[1,-2],    12,[2,3,6];
     -2,       4, [4,5];
     [-1,2,3,-4],2, [4,1]};

% This should be the polynomial-determinand of A
res_poly = polydet1(A)

function Y = polydet1(A)
    assert(ismatrix(A));
    s = size(A);
    assert( isequal(class(A),'cell') && (s(1) == s(2)) )
    
    if s(1) >= 10
        warning('Determinand of matrix with dimensions >9 requested. This can take a while. Are you shure you want to continue?');
        pause;
    end
    
    polymatrix = cell(size(A));
    for ii = 1:length(A(1,:))^2
        polymatrix{ii} = A{ii};
    end
    
    Y = polydet2(polymatrix);
end

function Y = polydet2(A)
    assert(iscell(A));
    s = size(A);
    assert( (s(1) == s(2)) || ~isequal(class(A{1}),'double') )
    
    m = s(1);
    
    if m == 2
        % 2x2 Matrix [A11 A12
        %             A21 A22]
        % return A11*A22 - A12*A21
        
        a = conv(A{1,1},A{2,2});
        b = conv(A{1,2},A{2,1});
        
        [a,b] = pad_left(a,b);
        Y = ( a - b );
    else
        Y = 0;
        for ii = 1:m
            B =     [A(1:ii-1,2:end);A(ii+1:end,2:end)];
            x =     (-1).^(ii-1) * conv(A{ii,1},polydet2(B));
            [x,Y] = pad_left(x,Y);
            Y =     Y + x;
        end
    end
end
