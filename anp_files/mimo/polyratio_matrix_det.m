function Y = polyratio_matrix_det(A)
    assert(ismatrix(A));
    s = size(A);
    assert( (s(1) == s(2)) || ~isequal(class(A),'polyratio') )
    
    if s(1) >= 10
        warning('Determinand of matrix with dimensions >9 requested. This can take a while. Are you shure you want to continue?');
        pause;
    end
    
    Y = polydet2(A);
    Y = Y.normalize;
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
        
        a = A{1,1}.mult(A{2,2});
        b = A{1,2}.mult(A{2,1});
        
        Y = a.sub(b);
    else
        Y = polyratio(0,1);
        for ii = 1:m
            B =     [A(1:ii-1,2:end);A(ii+1:end,2:end)];
            x =     polyratio((-1).^(ii-1),1).mult(A{ii,1}.mult(polydet2(B)));
            Y = Y.add(x);
        end
    end
end
