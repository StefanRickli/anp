s = tf('s');

tf11 = 1/(s+1);
tf12 = 2/(s+2);
tf13 = 3/(s+3);
tf21 = 4/(s+4);
tf22 = (s+1)/(s*(s-1));
tf23 = 6/(s-10);
tf31 = 7/(s+5);
tf32 = 8/(s+6);
tf33 = 9/(s+10);

G = [tf11,tf12,tf13;
     tf21,tf22,tf23;
     tf31,tf32,tf33];

res_poly = polydet1(G)

function Y = polydet1(A)
    assert(ismatrix(A));
    s = size(A);
    assert( (s(1) == s(2)) || ~isequal(class(A),'tf') )
    
    if s(1) >= 10
        warning('Determinand of matrix with dimensions >9 requested. This can take a while. Are you shure you want to continue?');
        pause;
    end
    
    polymatrix = cell(size(A));
    for ii = 1:length(A(1,:)).^2
        polymatrix{ii} = polyratio(A.Numerator{ii},A.Denominator{ii});
    end
    
    Y = polydet2(polymatrix);
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
