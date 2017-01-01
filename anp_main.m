function [] = anp_main(varargin)
    %close all;
    %clear classes;
    
    global debug debug_text;
    debug = false;
    debug_text = false;

    addpath('anp_files');
    addpath('anp_icons');

    args =                      anp_parse_arguments(varargin{:});

    % if argument checking has passed, get previous GUI and TF-Processing
    % objects back or re-initialize them if they have been deleted in the
    % meantime
    persistent h_tf_processor %h_gui

    h_gui = anp_gui;

    if isempty(h_tf_processor) || ~isvalid(h_tf_processor)
        h_tf_processor = anp_tf_processor;
    else
        fprintf('Reusing old TF processor instance\n');
    end

    args.gui_handle =           h_gui;
    args.processor_handle =     h_tf_processor;
    args.tf_poles =             roots(args.tf_obj.Numerator{1})';
    args.tf_zeros =             roots(args.tf_obj.Denominator{1})';

    h_tf_processor.init_params(args);

    h_gui.init_params(args);

    h_gui.init_visuals();

end

