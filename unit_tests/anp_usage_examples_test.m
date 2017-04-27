cd('..');
s = tf('s');

%% no arguments
out = anp_main('trigger_step','return_handle','cleanup_after_error');
drawnow;
delete(out{1});
delete(out{2});

%% radius specified
out = anp_main('trigger_step','return_handle','cleanup_after_error','R',30);
drawnow;
delete(out{1});
delete(out{2});

%% radius and 'right_dims' specified
out = anp_main('trigger_step','return_handle','cleanup_after_error','R',30,'right_dims',[.05 .05]);
drawnow;
delete(out{1});
delete(out{2});

%% zero and poles specified
out = anp_main([],[-1,-1],'trigger_step','return_handle','cleanup_after_error');
drawnow;
delete(out{1});
delete(out{2});

%% unstable system
unstable_system = tf(poly(-3),poly([1,-2]));
out = anp_main(unstable_system,'trigger_step','return_handle','cleanup_after_error');
drawnow;
delete(out{1});
delete(out{2});

%% unstable system, radius specified
unstable_system = tf(poly(-3),poly([1,-2]));
out = anp_main(unstable_system,'trigger_step','return_handle','cleanup_after_error','R',50);
drawnow;
delete(out{1});
delete(out{2});

%% large system
out = anp_main((10*(s+1)*(2*s+1))/((100*s+1)*(20*s+1)*(10*s+1)*(0.5*s+1)),'trigger_step','return_handle','cleanup_after_error');
drawnow;
delete(out{1});
delete(out{2});

%% large system, right_dims, right_x0 specified
out = anp_main((10*(s+1)*(2*s+1))/((100*s+1)*(20*s+1)*(10*s+1)*(0.5*s+1)),'trigger_step','return_handle','cleanup_after_error','right_dims',[1 1],'right_x0',[0 0]);
drawnow;
delete(out{1});
delete(out{2});

%% system with outliers 1
zeros = -30;
poles = [-2+1i,-2-1i];
out = anp_main(zeros,poles,'trigger_step','return_handle','cleanup_after_error');
drawnow;
delete(out{1});
delete(out{2});

%% system with outliers 2
zeros = -40;
poles = [-2+1i,-2-1i];
out = anp_main(zeros,poles,'trigger_step','return_handle','cleanup_after_error');
drawnow;
delete(out{1});
delete(out{2});
