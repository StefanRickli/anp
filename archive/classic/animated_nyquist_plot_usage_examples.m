% theory about Nyquist plot in Signals and Systems II, chapter 5
% download at http://control.ee.ethz.ch/~sigsys/ at the bottom

% usage (see animated_nyquist_plot.m file for full documentation):
% animated_nyquist_plot(tf)
% animated_nyquist_plot([zeros],[poles])
% animated_nyquist_plot(__,Name,Value)

fprintf('This script is not meant to run through, it might be tedious to abort!\n')
fprintf('Instead, execute the examples below by selecting the desired lines and hit ''F9''\n');
return;

% default is zeros = [-3] and poles = [-1+i,-1-i,-2]
% animation, contribution part
animated_nyquist_plot %#ok<UNRCH>

% change simple parameters like radius of D-curve and window dimensions of
% plots
animated_nyquist_plot('R',30);
animated_nyquist_plot('R',30,'right_dims',[.05 .05]);

% specify zeros and poles explicitly, SigSys II Ex from slide 5.31
animated_nyquist_plot([],[-1,-1])

% what happens in origin?
unstable_system = tf(poly([-3]),poly([1,-2]))
nyquist(unstable_system) % can't tell...

animated_nyquist_plot(unstable_system)
animated_nyquist_plot(unstable_system,'R',50)

% SigSys II Ex from slide 5.43
s = tf('s')
ex2 = (10*(s+1)*(2*s+1))/((100*s+1)*(20*s+1)*(10*s+1)*(0.5*s+1))
nyquist(ex2)
animated_nyquist_plot(ex2)
animated_nyquist_plot(ex2,'duration',60) % takes waay longer to compute but has better spatial resolution at the beginning
animated_nyquist_plot(ex2,'right_dims',[200 100],'right_x0',[-50 0])
animated_nyquist_plot(ex2,'right_dims',[20 10],'right_x0',[-5 0])
animated_nyquist_plot(ex2,'right_dims',[5 2],'right_x0',[-1 0])

% system with outliers, watch the zoom level in the left plot
zeros = [-30]
poles = [-2+i,-2-i]
animated_nyquist_plot(zeros,poles)

zeros = [-40]
poles = [-2+i,-2-i]
animated_nyquist_plot(zeros,poles)
