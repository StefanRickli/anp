
s = tf('s');

tf11 = 1/(s^2+4);
tf12 = 2/(s+2);
tf13 = 3/(s+3);
tf21 = 4/(s+4);
tf22 = (s+1)/(s-1);
tf23 = 6/(s-10);
tf31 = (s-5)/(s+5);
tf32 = 8/(s^2+2*s+6);
tf33 = 9/(s+10);

G = [tf11,tf12,tf13;
     tf21,tf22,tf23;
     tf31,tf32,tf33];

H = tf(ss([1,0.1;0.2,-2],[1,0.1;0.2,0.5],[1,0;0,1],[0,0;0,0]));

anp_main(H)

% polymatrix = cell(size(G));
% for ii = 1:length(G(1,:))^2
%     polymatrix{ii} = polyratio(G.Numerator{ii},G.Denominator{ii});
% end
% 
% for ii = 1:length(G(1,:))
%     for jj = 1:length(G(1,:))
%         disp(polymatrix{ii,jj})
%     end
% end
% 
% for ii = 1:length(G(1,:))
%     polymatrix{ii,ii} = polymatrix{ii,ii}.add(polyratio(1,1));
% end
% 
% for ii = 1:length(G(1,:))
%     for jj = 1:length(G(1,:))
%         disp(polymatrix{ii,jj})
%     end
% end
% 
% det_I_plus_L = polyratio_matrix_det(polymatrix);
% [z,p] = det_I_plus_L.reduce;
% 
% anp_main(tf(det_I_plus_L.num,det_I_plus_L.denom));