% theory about Nyquist plot in Signals and Systems II, chapter 5
% download at http://control.ee.ethz.ch/~sigsys/ at the bottom

% usage (see anp_main.m file for full documentation):
% anp_main(tf)
% anp_main([zeros],[poles])
% anp_main(__,Name,Value)

fprintf('This script is not meant to run through, it might be tedious to abort!\n')
fprintf('Instead, execute the examples below by selecting the desired lines and hit ''F9''\n');
return;

% default is zeros = [-3] and poles = [-1+i,-1-i,-2]
% animation, contribution part
anp_main %#ok<UNRCH>

% change simple parameters like radius of D-curve and window dimensions of
% plots
anp_main('R',30);
anp_main('R',30,'right_dims',[.05 .05]);
anp_main('R',30,'right_dims',[.001 .001]);

% specify zeros and poles explicitly, SigSys II Ex from slide 5.31
anp_main([],[-1,-1])

% what happens in origin?
unstable_system = tf(poly([-3]),poly([1,-2]))
nyquist(unstable_system) % can't tell...

anp_main(unstable_system)
anp_main(unstable_system,'R',50)

% SigSys II Ex from slide 5.43
s = tf('s')
ex2 = (10*(s+1)*(2*s+1))/((100*s+1)*(20*s+1)*(10*s+1)*(0.5*s+1))
nyquist(ex2)
anp_main(ex2)
anp_main(ex2,'right_dims',[1 1],'right_x0',[0 0])
anp_main(ex2,'right_dims',[0.01 0.01],'right_x0',[-0.001 0])
anp_main(ex2,'right_dims',[0.005 0.005],'right_x0',[-0.001 0])

% system with outliers, watch the zoom level in the left plot
zeros = [-30]
poles = [-2+i,-2-i]
anp_main(zeros,poles)

zeros = [-40]
poles = [-2+i,-2-i]
anp_main(zeros,poles)
