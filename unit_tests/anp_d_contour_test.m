cd('..');

%% no arguments
out = anp_main('trigger_step','return_handle','cleanup_after_error');
drawnow;
delete(out{1});
delete(out{2});

%% p pos overlapping
out = anp_main([],[0.2i],'trigger_step','return_handle','cleanup_after_error');
drawnow;
delete(out{1});
delete(out{2});

%% p neg overlapping
out = anp_main([],[-0.2i],'trigger_step','return_handle','cleanup_after_error');
drawnow;
delete(out{1});
delete(out{2});

%% p pos on origin
out = anp_main([],[0.25i],'trigger_step','return_handle','cleanup_after_error');
drawnow;
delete(out{1});
delete(out{2});

%% p neg on origin
out = anp_main([],[-0.25i],'trigger_step','return_handle','cleanup_after_error');
drawnow;
delete(out{1});
delete(out{2});

%% z pos overlapping
out = anp_main([0.1i],[-2],'trigger_step','return_handle','cleanup_after_error');
drawnow;
delete(out{1});
delete(out{2});

%% z neg overlapping
out = anp_main([-0.1i],[-2],'trigger_step','return_handle','cleanup_after_error');
drawnow;
delete(out{1});
delete(out{2});

%% z pos on origin
out = anp_main([0.125],[-2],'trigger_step','return_handle','cleanup_after_error');
drawnow;
delete(out{1});
delete(out{2});

%% z neg on origin
out = anp_main([-0.125],[-2],'trigger_step','return_handle','cleanup_after_error');
drawnow;
delete(out{1});
delete(out{2});

%% p neg on origin + p pos
out = anp_main([],[1i,-0.25i],'trigger_step','return_handle','cleanup_after_error');
drawnow;
delete(out{1});
delete(out{2});

%% p pos on origin + p pos
out = anp_main([],[0.25i,2i],'trigger_step','return_handle','cleanup_after_error');
drawnow;
delete(out{1});
delete(out{2});

%% p pos on origin + p neg
out = anp_main([],[0.25i,-2i],'trigger_step','return_handle','cleanup_after_error');
drawnow;
delete(out{1});
delete(out{2});

%% p pos overlapping + p pos
out = anp_main([],[0.2i,2i],'trigger_step','return_handle','cleanup_after_error');
drawnow;
delete(out{1});
delete(out{2});

%% p pos overlapping + p neg
out = anp_main([],[0.2i,-2i],'trigger_step','return_handle','cleanup_after_error');
drawnow;
delete(out{1});
delete(out{2});

%% p pos overlapping + z pos
out = anp_main([1i],[0.2i],'trigger_step','return_handle','cleanup_after_error');
drawnow;
delete(out{1});
delete(out{2});

%% p pos overlapping + z neg
out = anp_main([-1i],[0.2i],'trigger_step','return_handle','cleanup_after_error');
drawnow;
delete(out{1});
delete(out{2});

%% p neg overlapping + z pos
out = anp_main([1i],[-0.2i],'trigger_step','return_handle','cleanup_after_error');
drawnow;
delete(out{1});
delete(out{2});

%% p neg overlapping + z neg
out = anp_main([-1i],[-0.2i],'trigger_step','return_handle','cleanup_after_error');
drawnow;
delete(out{1});
delete(out{2});

%% p pos on origin + z pos
out = anp_main([1i],[0.25i],'trigger_step','return_handle','cleanup_after_error');
drawnow;
delete(out{1});
delete(out{2});

%% p pos on origin + z neg
out = anp_main([-1i],[0.25i],'trigger_step','return_handle','cleanup_after_error');
drawnow;
delete(out{1});
delete(out{2});

%% p neg on origin + z pos
out = anp_main([1i],[-0.25i],'trigger_step','return_handle','cleanup_after_error');
drawnow;
delete(out{1});
delete(out{2});

%% p neg on origin + z neg
out = anp_main([-1i],[-0.25i],'trigger_step','return_handle','cleanup_after_error');
drawnow;
delete(out{1});
delete(out{2});

%% z pos on origin + p pos
out = anp_main([0.125i],[1i],'trigger_step','return_handle','cleanup_after_error');
drawnow;
delete(out{1});
delete(out{2});

%% z pos on origin + p neg
out = anp_main([0.125i],[-1i],'trigger_step','return_handle','cleanup_after_error');
drawnow;
delete(out{1});
delete(out{2});

%% z neg on origin + p pos
out = anp_main([-0.125i],[1i],'trigger_step','return_handle','cleanup_after_error');
drawnow;
delete(out{1});
delete(out{2});

%% z neg on origin + p neg
out = anp_main([-0.125i],[-1i],'trigger_step','return_handle','cleanup_after_error');
drawnow;
delete(out{1});
delete(out{2});

%% z pos overlapping + p pos
out = anp_main([0.1i],[1i],'trigger_step','return_handle','cleanup_after_error');
drawnow;
delete(out{1});
delete(out{2});

%% z pos overlapping + p neg
out = anp_main([0.1i],[-1i],'trigger_step','return_handle','cleanup_after_error');
drawnow;
delete(out{1});
delete(out{2});

%% z neg overlapping + p pos
out = anp_main([-0.1i],[1i],'trigger_step','return_handle','cleanup_after_error');
drawnow;
delete(out{1});
delete(out{2});

%% z neg overlapping + p neg
out = anp_main([-0.1i],[-1i],'trigger_step','return_handle','cleanup_after_error');
drawnow;
delete(out{1});
delete(out{2});

%% z pos overlapping + z pos
out = anp_main([0.1i,2i],[-2,-4],'trigger_step','return_handle','cleanup_after_error');
drawnow;
delete(out{1});
delete(out{2});

%% z pos overlapping + z neg
out = anp_main([0.1i,-2i],[-2,-4],'trigger_step','return_handle','cleanup_after_error');
drawnow;
delete(out{1});
delete(out{2});

%% z neg overlapping + z pos
out = anp_main([-0.1i,2i],[-2,-4],'trigger_step','return_handle','cleanup_after_error');
drawnow;
delete(out{1});
delete(out{2});

%% z neg overlapping + z neg
out = anp_main([-0.1i,-2i],[-2,-4],'trigger_step','return_handle','cleanup_after_error');
drawnow;
delete(out{1});
delete(out{2});

%% z pos on origin + z pos
out = anp_main([0.125i,2i],[-2,-4],'trigger_step','return_handle','cleanup_after_error');
drawnow;
delete(out{1});
delete(out{2});

%% z pos on origin + z neg
out = anp_main([0.125i,-2i],[-2,-4],'trigger_step','return_handle','cleanup_after_error');
drawnow;
delete(out{1});
delete(out{2});

%% z neg on origin + z pos
out = anp_main([-0.125i,2i],[-2,-4],'trigger_step','return_handle','cleanup_after_error');
drawnow;
delete(out{1});
delete(out{2});

%% z neg on origin + z neg
out = anp_main([0.125i,2i],[-2,-4],'trigger_step','return_handle','cleanup_after_error');
drawnow;
delete(out{1});
delete(out{2});

