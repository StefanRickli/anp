function [] = anp_main(varargin)
%     try
        global debug;
        debug = true;

        addpath('anp_files');
        addpath('anp_icons');

        checked_args = anp_parse_arguments(varargin{:});

        % if argument checking has passed, get previous GUI and TF-Processing
        % objects back or re-initialize them if they have been deleted in the
        % meantime
        persistent h_gui h_tf_processor

        if isempty(h_gui) || ~isvalid(h_gui)
            h_gui = anp_gui;
        else
            fprintf('Reusing old GUI instance\n');
        end
        if isempty(h_tf_processor) || ~isvalid(h_tf_processor)
            h_tf_processor = anp_tf_processor;
        else
            fprintf('Reusing old TF processor instance\n');
        end
            
        
        processor_args.tf_obj =     checked_args.tf_obj;
        processor_args.radii =      checked_args.radii;
        processor_args.angles =     checked_args.angles;
        processor_args.weights =    checked_args.weights;
        
        h_tf_processor.set_all_params(processor_args);
        
        gui_args = struct('z_x0',[],...
                          'z_dims',[],...
                          'z_auto_lims',[],...
                          'w_x0',[],...
                          'w_dims',[],...
                          'w_auto_lims',[],...
                          'poles',[],...
                          'zeros',[],...
                          'radii',[],...
                          'border',[],...
                          'plot_size',[],...
                          'trail_length',[]);
        
        h_gui.set_anp_tf_processor(h_tf_processor);
        gui_args.z_x0 =         checked_args.z_plot_x0;
        gui_args.z_dims =       checked_args.z_plot_dims;
        gui_args.z_auto_lims =  checked_args.z_plot_auto_lims;
        gui_args.w_x0 =         checked_args.w_plot_x0;
        gui_args.w_dims =       checked_args.w_plot_dims;
        gui_args.w_auto_lims =  checked_args.w_plot_auto_lims;
        gui_args.poles =        roots(checked_args.tf_obj.Numerator{1})';
        gui_args.zeros =        roots(checked_args.tf_obj.Denominator{1})';
        gui_args.radii.R =        5;% TODO checked_args.radii;
        gui_args.border =       checked_args.border;
        gui_args.plot_size =    checked_args.plot_size;
        gui_args.trail_length =    checked_args.trail_length;
        
        h_gui.set_z_plot_limits(gui_args);
        h_gui.set_w_plot_limits(gui_args);
        h_gui.set_poles_zeros(gui_args);
        % TODO radii comes from tf_processor!
        h_gui.set_radii(gui_args);
        h_gui.set_window_props(gui_args);
        h_gui.set_plot_props(gui_args);
        
        h_gui.init_visuals();
        
        
        
        
        %rmpath('anp_files');
%     catch err
%         disp(err.message);
%         keyboard;
%         clear all;
%         close all;
%     end
end

