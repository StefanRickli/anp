% For the theory behind Nyquist plots and their usage refer to
% http://lpsa.swarthmore.edu/Nyquist/Nyquist.html
% The program "NyquistGUI" there is also a nice alternative to my program.
% 
% If you rather like to look the theory up in a book, then I can recommend
% "Feedback Control of Dynamic Systems" by Franklin, Powell, Emami-Naeini

% Usage (see anp_main.m file for full documentation of parameters):
% anp_main(tf)
% anp_main([zeros],[poles])
% anp_main(__,Name,Value)

fprintf('This script is not meant to run through, it might be tedious to abort!\n')
fprintf('Instead, execute the examples below by selecting the desired lines and hit ''F9''\n');
return;

% Absence of arguments:
% default to zeros = [-0.7] and poles = [-1+i,-1-i,-2,-3]
anp_main                                                        %#ok<UNRCH>

% change simple parameters like radius of D-curve and window dimensions of
% plots
anp_main('R',30);
anp_main('R',30,'right_dims',[.05 .05]);
anp_main('R',30,'right_dims',[.001 .001]);

% specify zeros and poles explicitly
anp_main([],[-1,-1])

% what happens in origin?
unstable_system = tf(poly([-3]),poly([1,-2]))
nyquist(unstable_system) % can't tell...

anp_main(unstable_system)
anp_main(unstable_system,'R',50)

% Large system
s = tf('s')
sys = (10*(s+1)*(2*s+1))/((100*s+1)*(20*s+1)*(10*s+1)*(0.5*s+1))
nyquist(sys)
anp_main(sys)
anp_main(sys,'right_dims',[1 1],'right_x0',[0 0])
anp_main(sys,'right_dims',[0.01 0.01],'right_x0',[-0.001 0])
anp_main(sys,'right_dims',[0.005 0.005],'right_x0',[-0.001 0])

% system with outliers, watch the zoom level in the left plot once the zero
% is far away on the left
zeros = [-30]
poles = [-2+i,-2-i]
anp_main(zeros,poles)

zeros = [-40]
poles = [-2+i,-2-i]
anp_main(zeros,poles)

% random SISO 4-state system
G = tf(rss(4,1,1));
anp_main(G);

% random 2x2 MIMO (3 states per tf) system, the plot shows det(I - G(s))
G = tf(rss(3,2,2));
anp_main(G);

% random 4x4 MIMO (2 state per tf) system, the plot shows det(I - G(s))
G = tf(rss(2,4,4));
anp_main(G);
