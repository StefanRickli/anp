s = tf('s');

%% 3x3 system with lots of different integer poles
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

out = anp_main(G,'trigger_step','return_handle','cleanup_after_error');
drawnow;
delete(out{1});
delete(out{2});

%% 3x3 system out of state space system
G = tf(ss([-2,1,0.1;0.2,-2,-0.1;0.05,-0.5,-1],[1,0.1;0.2,0.5;0.1,-0.7],[1,0,0;0,1,0],[0,0;0,0]));

out = anp_main(G,'trigger_step','return_handle','cleanup_after_error');
drawnow;
delete(out{1});
delete(out{2});
