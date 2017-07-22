cd('..');

%% random SISO 1-state system
G = tf(rss(1,1,1));

out = anp_main(G,'trigger_step','return_handle','cleanup_after_error');
drawnow;
delete(out{1});
delete(out{2});

%% random SISO 2-state system
G = tf(rss(2,1,1));

out = anp_main(G,'trigger_step','return_handle','cleanup_after_error');
drawnow;
delete(out{1});
delete(out{2});

%% random SISO 3-state system
G = tf(rss(3,1,1));

out = anp_main(G,'trigger_step','return_handle','cleanup_after_error');
drawnow;
delete(out{1});
delete(out{2});

%% random SISO 4-state system
G = tf(rss(4,1,1));

out = anp_main(G,'trigger_step','return_handle','cleanup_after_error');
drawnow;
delete(out{1});
delete(out{2});

%% random 2x2 MIMO (1 state per tf) system
G = tf(rss(1,2,2));

out = anp_main(G,'trigger_step','return_handle','cleanup_after_error');
drawnow;
delete(out{1});
delete(out{2});

%% random 2x2 MIMO (2 states per tf) system
G = tf(rss(2,2,2));

out = anp_main(G,'trigger_step','return_handle','cleanup_after_error');
drawnow;
delete(out{1});
delete(out{2});

%% random 2x2 MIMO (3 states per tf) system
G = tf(rss(3,2,2));

out = anp_main(G,'trigger_step','return_handle','cleanup_after_error');
drawnow;
delete(out{1});
delete(out{2});

%% random 3x3 MIMO (1 state per tf) system
G = tf(rss(1,3,3));

out = anp_main(G,'trigger_step','return_handle','cleanup_after_error');
drawnow;
delete(out{1});
delete(out{2});

%% random 3x3 MIMO (2 states per tf) system
G = tf(rss(2,3,3));

out = anp_main(G,'trigger_step','return_handle','cleanup_after_error');
drawnow;
delete(out{1});
delete(out{2});

%% random 4x4 MIMO (1 state per tf) system
G = tf(rss(1,4,4));

out = anp_main(G,'trigger_step','return_handle','cleanup_after_error');
drawnow;
delete(out{1});
delete(out{2});

%% random 4x4 MIMO (2 state per tf) system
G = tf(rss(2,4,4));

out = anp_main(G,'trigger_step','return_handle','cleanup_after_error');
drawnow;
delete(out{1});
delete(out{2});
