%   class: anp_gui
%   --------------
%   
%   This class is all about displaying the data that a tf_processor object
%   previously calculated.
%   
%   On creation it instantiates a figure containing a left and right
%   subplot and some space below to put text about the transfer function
%   that is going to be plotted.
%   
%   After that it takes the data from a tf_processor object, calculates
%   display and animation parameters and draws it to the preallocated
%   graphs and text boxes.
%
%   The figure contains a custom toolbar with buttons that control the
%   stepping of the animation.
%
%   paper reference: p.41
classdef anp_gui < handle
    
    % ---------------------------------------------------------------------
    % PRIVATE Variables
    % ---------------------------------------------------------------------
        
    properties(SetAccess = private)
        % g: general properties
        
        g_uid               % Java String,      uniqe identifier for every object instance
        
        
        % s: state
        
        s_gui_nominal       logical             % do we have everything together before we can draw?
        s_draw_allowed      logical             % set this to false before changing data or parameters, it forces the draw-function to stop once it's done with an iteration
        s_draw_busy         logical             % we're updating the plot at the moment, parameter or data change not allowed
        s_check_limits      logical             % prompt the draw-function to update axis limit related data before it makes the next iteration
        
        
        % h: handles
        
        h_anp_tf_processor                      % remember the handle to the calculation object
        
        h_fig               matlab.ui.Figure                            % main figure handle
        h_sub1              matlab.graphics.axis.Axes                   % left (z-)plot
        h_sub2              matlab.graphics.axis.Axes                   % right (w-)plot
        h_zoom              matlab.graphics.interaction.internal.zoom   % oject needed for zooming related callbacks
        h_pan               matlab.graphics.interaction.internal.pan    % oject needed for panning related callbacks
        
        h_text_p_annot      matlab.graphics.shape.TextBox           % handles to the textboxes of the pole-contributions
        h_text_z_annot      matlab.graphics.shape.TextBox           % handles to the textboxes of the zero-contributions
        h_text_res_annot    matlab.graphics.shape.TextBox           % handles to the textboxes of the result
        
        h_z_plot_full       matlab.graphics.chart.primitive.Line    % left plot full D-shape
        h_z_plot_trail      matlab.graphics.chart.primitive.Line    % left plot yellow/orange trail of arrow
        h_z_plot_arrow      matlab.graphics.shape.Arrow             % left plot arrow
        h_z_pz_objcts       % annotations of type crosses, circles and arrows with text
        
        h_w_plot_full       matlab.graphics.chart.primitive.Line
        h_w_plot_trail      matlab.graphics.chart.primitive.Line
        h_w_plot_arrow      matlab.graphics.shape.Arrow
        
        
        % ui: GUI elements
        
        ui_icons            double                          % 4D array (x,y,color,index) that holds the custom icons
        ui_toolbar          matlab.ui.container.Toolbar
        ui_run_switches     matlab.ui.container.toolbar.ToggleTool
        ui_sw_pause         matlab.ui.container.toolbar.ToggleTool
        ui_btn_prev         matlab.ui.container.toolbar.PushTool
        ui_btn_next         matlab.ui.container.toolbar.PushTool
        
        
        % a: animation
        
        a_direction         double          % determines movement speed of the arrow. a_direction = {-3,-1,1,3}
        a_time_ii           double          % remember the position of the arrow

        
        % w: window properties
        
        w_fig_position      double          % vector containing the figure window position on the display
        w_sub1_position     double          % vector conataining the left plot's position withing the figure
        w_sub2_position     double          % vector conataining the right plot's position withing the figure
        
        w_border                    double  % how much border to the left, right and top from the plots? [pixel]
        w_plot_size                 double  % [pixel]
        w_plot_width_frac           double  % [1]
        w_border_horizontal_frac    double  % [1]
        w_border_vertical_frac      double  % [1]
        w_annotation_start_frac     double  % where (vertically) on the figure do the annotations start? [1]
        w_annotation_textbox_frac   double  % vertical separation between annotations [1]
        
        
        % p: plot properties (z: left, w=f(z): right)
        
        p_n_time_steps      double          % how many steps in the animation? limit for a_time_ii
        p_oversampling_factor double        % how many more data than animation steps are there?
        p_n_data_points       double        % basically time steps * oversampling factor
        p_trail_length      double          % which percentage of data points make the trail?
        p_n_trail           double          % of how many data points is the trail comprised?
        
        p_z_x0              double          % left plot: requested center location
        p_z_dims            double          % left plot: requested height and width
        p_z_auto_lims       logical         % left plot: ignore p_z_x0 and p_z_dims, calculate xlim and ylim automatically
        p_z_xlim            double          % left plot: actual xlim
        p_z_xspan           double          % left plot: data range in x-direction
        p_z_ylim            double          % left plot: actual ylim
        p_z_yspan           double          % left plot: data range in y-direction
        p_z_arrow_length    double          % left plot: length of annotation as fraction of the figure [1]
        
        p_w_x0              double
        p_w_dims            double
        p_w_auto_lims       logical
        p_w_xlim            double
        p_w_xspan           double
        p_w_ylim            double
        p_w_yspan           double
        p_w_arrow_length    double
        
        
        % d: data
        
        d_n_poles               double      % how many poles does the transfer function have?
        d_n_zeros               double      % how many zeros does the transfer function have?
        d_poles                 double      % array holding the pole locations
        d_zeros                 double      % array holding the zero locations
        d_R                     double      % radius of the half-circle
        
        d_t_values              double      % UNUSED, time steps
        d_t_oversampled         double      % UNUSED
        d_t_trails                          % cell array holding the index arrays that address the different trail elements at a certain animation step
        
        d_z_values              double      % array that holds the raw data for the left plot
        d_z_values_truncated    double      % array, holds the displayed data of the left plot. it might be altered in order to fit the x- and ylim.
        
        d_w_values              double
        d_w_values_truncated    double        
    end
    
    
    % ---------------------------------------------------------------------
    % PUBLIC methods
    % ---------------------------------------------------------------------
    
    methods(Access = public)
        
        function this = anp_gui()
            % Instantiates a figure window with two subplots. Initializes state variables and the custom toolbar.
            
            % set a unique ID for every instance of the GUI and find out
            % which Matlab version is running.
            this.g_uid =        java.util.UUID.randomUUID.toString;
            
            % Initialize the state variables
            this.s_gui_nominal = false;
            this.s_draw_allowed = false;
            this.s_draw_busy =  false;
            this.a_time_ii =    1;
            this.a_direction =  1;
            
            % Instantiate the figure and its subplots
            this.h_fig =        figure('IntegerHandle','off','CloseRequestFcn',@this.on_figure_delete);
            this.h_sub1 =       subplot(1,2,1);
            hold on;
            this.h_sub2 =       subplot(1,2,2);
            hold on;
            this.h_fig.HandleVisibility = 'off';
            
            % Initialize the custom toolbar
            this.ui_toolbar =   findall(this.h_fig,'Type','uitoolbar');
            this.load_icons();
            this.ui_run_switches(1) =   uitoggletool(this.ui_toolbar,'CData',this.ui_icons(:,:,:,1),'TooltipString','Fast Backward','OnCallback',{@this.cb_run,-3},'OffCallback',{@this.cb_run_switch_off,-3},'Enable','off','Separator','on');
            this.ui_run_switches(2) =   uitoggletool(this.ui_toolbar,'CData',this.ui_icons(:,:,:,2),'TooltipString','Play Backward','OnCallback',{@this.cb_run,-1},'OffCallback',{@this.cb_run_switch_off,-1},'Enable','off');
            this.ui_btn_prev =          uipushtool  (this.ui_toolbar,'CData',this.ui_icons(:,:,:,3),'TooltipString','Step Back',    'ClickedCallback',{@this.cb_step,-1},'Enable','off');
            this.ui_sw_pause =          uitoggletool(this.ui_toolbar,'CData',this.ui_icons(:,:,:,4),'TooltipString','Pause',        'ClickedCallback',{@this.cb_pause},'Enable','off','State','on');
            this.ui_btn_next =          uipushtool  (this.ui_toolbar,'CData',this.ui_icons(:,:,:,5),'TooltipString','Step Forward', 'ClickedCallback',{@this.cb_step,1},'Enable','off');
            this.ui_run_switches(3) =   uitoggletool(this.ui_toolbar,'CData',this.ui_icons(:,:,:,6),'TooltipString','Play Forward', 'OnCallback',{@this.cb_run,1},'OffCallback',{@this.cb_run_switch_off,1},'Enable','off');
            this.ui_run_switches(4) =   uitoggletool(this.ui_toolbar,'CData',this.ui_icons(:,:,:,7),'TooltipString','Fast Forward', 'OnCallback',{@this.cb_run,3},'OffCallback',{@this.cb_run_switch_off,3},'Enable','off');
            
            tools.dbg('anp_gui[constructor]:\t%s: Instance created.\n',this.g_uid);
        end
        
        function [] = init_params(this,new_props)
            % Populate the variables that hold the relevant properties and
            % data for the GUI
            
            % TODO  if we're reusing this instance, check for dirty values
            %       and later redraw
            
            this.set_z_plot_limits(new_props);
            this.set_w_plot_limits(new_props);
            this.set_poles_zeros(new_props);
            this.set_window_props(new_props);
            this.set_plot_props(new_props);
            
            this.set_anp_tf_processor(new_props);
        end
        
        function [] = set_window_props(this,new_props)
            % Populates the variables holding the main window properties with new values.
            
            this.w_border =     new_props.border;
            this.w_plot_size =  new_props.plot_size;
        end
        
        function [] = set_plot_props(this,new_props)
            % Populates the variables holding the general plot properties with new values.

            this.p_trail_length =     new_props.trail_length;
        end
        
        function [] = set_z_plot_limits(this,new_props)
            % Populate the axis parameters of the left plot with new values.
            
            % TODO check for proper array dimensions
            
            this.p_z_x0 =           new_props.z_plot_x0;
            this.p_z_dims =         new_props.z_plot_dims;
            this.p_z_auto_lims =    new_props.z_plot_auto_lims;
        end
        
        function [] = set_w_plot_limits(this,new_props)
            % Populate the axis parameters of the right plot with new values.
            
            % TODO check for proper array dimensions
            
            this.p_w_x0 =           new_props.w_plot_x0;
            this.p_w_dims =         new_props.w_plot_dims;
            this.p_w_auto_lims =    new_props.w_plot_auto_lims;
        end
        
        function [] = set_radii(this,new_props)
            % Populate variables holding geometric information with new values.
            
            % TODO check for proper array dimensions
            % TODO radii comes from tf_processor!
            this.d_radii =          new_props.radii;
        end
        
        function [] = set_poles_zeros(this,new_props)
            % Populates the variables containing information about the poles and zeros with new values.
            
            % TODO check for proper array dimensions
            
            this.d_poles =          new_props.tf_poles;
            this.d_zeros =          new_props.tf_zeros;
            
            this.d_n_poles =        length(new_props.tf_poles);
            this.d_n_zeros =        length(new_props.tf_zeros);
        end

        function [] = set_anp_tf_processor(this,new_props)
            % Changes the handle reference to a new tf_processor object.
            
            this.h_anp_tf_processor = new_props.processor_handle;
        end
        
        function [] = init_visuals(this)
            % Calculates everything needed for a first draw and performs it.
            
            % Get data and relevant display parameters from tf_processor.
            this.fetch_R();
            this.fetch_data();
            
            % Decide on layout parameters.
            this.calc_gui_positions();
            this.calc_axis_limits();
            this.calc_z_arrow_length();
            this.calc_w_arrow_length();
            this.calc_trail_indexes();
            
            % Calculate the data that is to be plotted. Compared to the raw
            % data from the tf_processor object, its data points could be
            % repositioned to stay within the plot.
            this.calc_truncated_z_values();
            this.calc_truncated_w_values();
            
            % Instatiate the rest of the objects in the figure and/or
            % configure them.
            this.draw_init_gui_statics();
            this.draw_init_gui_text_objects();
            this.draw_init_plot_axes();
            this.draw_init_line_plots();
            this.draw_init_plot_arrows();
            this.draw_init_z_plot_poles_zeros();
            
            % Since draw_one_frame only updates the trail arrow and the
            % trail itself, we must explicitly update the plot of the whole
            % functions.
            this.draw_update_full_z_plot();
            this.draw_update_full_w_plot();
            
            % Update trail, trail arrow and text boxes.
            this.draw_one_frame();
            
            % s_gui_nominal = true means that everything has been
            % initialized and the GUI can enter a nominal state.
            this.s_gui_nominal = true;
            
            % Until now the buttons were disabled. Make them clickable now.
            this.ui_control_enable();
        end
        
        
        function delete(this)
            % Destructor for the object. Closes the figure if it hasn't been already.
            
            tools.dbg('anp_gui[delete]\t%s: Deletion requested.\n',this.g_uid);
            
            if isvalid(this.h_fig)
                close(this.h_fig);
            end
        end
    end
    
    % ---------------------------------------------------------------------
    % PRIVATE methods
    % ---------------------------------------------------------------------
    
    methods(Access = private)
        function on_figure_delete(this,src,~) % ignored argument is 'evt'
            % Callback that issues the deletion of the anp_gui object that is holding the closed and now invalid figure window.
            
            tools.dbg('anp_gui[on_figure_delete]:\tDeleting figure window.\n');
            
            delete(src);
            this.delete();
        end
        
        function [] = load_icons(this)
            % Loads the PNG files of the custom toolbar icons into memory.
            
            this.ui_icons =             zeros(16,16,3,5);
            this.ui_icons(:,:,:,1) =    imread('Fast Backward.png');
            this.ui_icons(:,:,:,2) =    imread('Backward.png');
            this.ui_icons(:,:,:,3) =    imread('Step Backward.png');
            this.ui_icons(:,:,:,4) =    imread('Pause.png');
            this.ui_icons(:,:,:,5) =    imread('Step Forward.png');
            this.ui_icons(:,:,:,6) =    imread('Forward.png');
            this.ui_icons(:,:,:,7) =    imread('Fast Forward.png');
            this.ui_icons =             this.ui_icons/255;
        end
        
        function [] = fetch_R(this)
            % Gets the current value of the half-circle from the tf_processor object.
            
            assert(isvalid(this.h_anp_tf_processor));
            
            this.d_R =                  this.h_anp_tf_processor.get_R();
        end
        
        function [] = fetch_data(this)
            % Gets all data that could be relevant to the GUI from the tf_processor object.
            
            assert(isvalid(this.h_anp_tf_processor));
            
            time_data =                     this.h_anp_tf_processor.get_time_points();
            this.d_t_values =               time_data.time_points;
            this.d_t_oversampled =          time_data.data_points;
            this.p_n_time_steps =           time_data.time_props.n_time_steps;
            this.p_oversampling_factor =    time_data.time_props.oversampling_factor;
            this.p_n_data_points =          time_data.time_props.n_data_points;
            
            function_data =         this.h_anp_tf_processor.get_data();
            this.d_z_values =       function_data.z_values;
            this.d_w_values =       function_data.w_values;
        end
        
        % ---------------
        % CALC-Methods
        % ---------------
        
        function [] = calc_gui_positions(this)
            % Determines where within the figure the different elements should be.
            % The function does this first in units of pixel (for better
            % human readable code) and then converts the findings to
            % fractions like the .Position function expects.
            
            % Calculate the correct positions (relative to interior of the
            % figure).
            fig_plot_height =               2*this.w_border + this.w_plot_size; % plot + border above and below, [pixel]
            fig_annotation_textbox_height = 14;                                 % one box, [pixel]
            fig_annotation_height_sum =     (max(this.d_n_zeros,this.d_n_poles)+2) ...
                                              * fig_annotation_textbox_height;  % cumulative, [pixel]
            
            fig_height =            fig_plot_height ...
                                    + fig_annotation_height_sum;                % [pixel]
            fig_width =             3*this.w_border + 2*this.w_plot_size;       % [pixel]
            
            % fig.Position expects pixels as unit.
            this.w_fig_position =   [100, ...
                                     100, ...
                                     fig_width, ...
                                     fig_height];
            
            % As subplot.Postion expects fractions of the inside of the
            % figure, we recalculate the pixel values into fractions.
            this.w_border_horizontal_frac = this.w_border / fig_width;                                      % [1]
            this.w_border_vertical_frac =   this.w_border / fig_height;                                     % [1]
            this.w_plot_width_frac =        (fig_width - 3*this.w_border) / (2*fig_width);                  % [1]
            fig_plot_height_frac =          1 - (2*this.w_border + fig_annotation_height_sum) / fig_height; % [1]
            
            % subX.Position expects fractions of the inside of the figure
            % as unit.
            this.w_sub1_position =  [this.w_border_horizontal_frac, ...
                                     (1 - this.w_border_vertical_frac - fig_plot_height_frac), ...
                                     this.w_plot_width_frac, ...
                                     fig_plot_height_frac];
            this.w_sub2_position =  [(2*this.w_border_horizontal_frac + this.w_plot_width_frac), ...
                                     (1 - this.w_border_vertical_frac - fig_plot_height_frac), ...
                                     this.w_plot_width_frac, ...
                                     fig_plot_height_frac];
            
            % Text box object parameters:
            % Again we need to know relative (fractions) positions inside
            % the figure.
            this.w_annotation_start_frac =      (fig_annotation_height_sum + 10) / fig_height;  % [1]
            this.w_annotation_textbox_frac =    fig_annotation_textbox_height / fig_height;     % [1]
        end
        
        function [] = calc_axis_limits(this)
            % Sets the correct x- and ylims for the left and right plot.
            % The method respects the auto-flags and only calculates the
            % parameters if it's told to do so.
            % Otherwise it uses the user given limits.
            
            % TODO make sure that all props and data are set before calling
            % this method.
            
            assert(~any([isempty([this.d_zeros,this.d_poles]),...
                         isempty(this.d_R) || isnan(this.d_R),...
                         isempty(this.d_w_values),...
                         isempty(this.p_z_auto_lims),...
                         isempty(this.p_w_auto_lims)]));
            
            if this.p_z_auto_lims
                % try to find optimal axis limits for the left plot, or let
                % the user decide by setting either 'left_x0' or
                % 'left_dims' to a manual value
                
                [this.p_z_xlim,this.p_z_ylim] = anp_plot_auto_zoom_z([this.d_zeros,this.d_poles],this.d_R); % separate source file
            else
               this.p_z_xlim = [(this.p_z_x0(1) - this.p_z_dims(1)/2), ...
                                (this.p_z_x0(1) + this.p_z_dims(1)/2)];
               
               this.p_z_ylim = [(this.p_z_x0(2) - this.p_z_dims(2)/2), ...
                                (this.p_z_x0(2) + this.p_z_dims(2)/2)];
            end
            
            % Override too small axis limits by replacing them with the
            % span of the data.
            [this.p_z_xspan,this.p_z_yspan] = anp_plot_find_span([this.d_zeros,this.d_poles,this.d_z_values]); % separate source file
            if diff(this.p_z_xlim) < 100*eps
                this.p_z_xlim = this.p_z_xspan;
            end
            if diff(this.p_z_ylim) < 100*eps
                this.p_z_ylim = this.p_z_yspan;
            end
            
            
            if this.p_w_auto_lims
                % try to find optimal axis limits for the right plot, or
                % let the user decide by setting either 'right_x0' or
                % 'right_dims' to a manual value
                
                [this.p_w_xlim,this.p_w_ylim] = anp_plot_auto_zoom_w(this.d_w_values); % separate source file
            else
               this.p_w_xlim = [(this.p_w_x0(1) - this.p_w_dims(1)/2), ...
                                (this.p_w_x0(1) + this.p_w_dims(1)/2)];
               
               this.p_w_ylim = [(this.p_w_x0(2) - this.p_w_dims(2)/2), ...
                                (this.p_w_x0(2) + this.p_w_dims(2)/2)];
            end
            
            % Override too small axis limits by replacing them with the
            % span of the data.
            [this.p_w_xspan,this.p_w_yspan] = anp_plot_find_span(this.d_w_values); % separate source file
            if diff(this.p_w_xlim) < 100*eps
                this.p_w_xlim = this.p_w_xspan;
            end
            if diff(this.p_w_xlim) < 100*eps
                this.p_w_ylim = this.p_w_yspan;
            end
        end
        
        function [] = calc_z_arrow_length(this)
            % Calculates how long the trail arrow should be based on the current axis limits.
            
            in_axis_width =         diff(this.p_z_xlim);                            % [1]
            in_axis_height =        diff(this.p_z_ylim);                            % [1]
            this.p_z_arrow_length = 0.04 * sqrt(in_axis_width^2 + in_axis_height^2);% [1]
            tools.dbg('anp_gui[calc_z_arrow_length]:\t%.5f\n',this.p_z_arrow_length);
        end
        
        function [] = calc_w_arrow_length(this)
            % Calculates how long the trail arrow should be based on the current axis limits.
            
            out_axis_width =        diff(this.p_w_xlim);                                % [1]
            out_axis_height =       diff(this.p_w_ylim);                                % [1]
            this.p_w_arrow_length = 0.04 * sqrt(out_axis_width^2 + out_axis_height^2);  % [1]
            tools.dbg('anp_gui[calc_w_arrow_length]:\t%.5f\n',this.p_w_arrow_length);
        end
        
        function [] = calc_truncated_z_values(this)
            % Wrapper function for this.trunc(). Moves values that are outside the coordinate system into the graph so we get a continuous closed curve.
            
            this.d_z_values_truncated =  this.trunc(this.d_z_values,this.p_z_xlim,this.p_z_ylim);
        end
        
        function [] = calc_truncated_w_values(this)
            % Wrapper function for this.trunc(). Moves values that are outside the coordinate system into the graph so we get a continuous closed curve.
            
            this.d_w_values_truncated =  this.trunc(this.d_w_values,this.p_w_xlim,this.p_w_ylim);
        end
        
        function [] = calc_trail_indexes(this)
            % Prepare t-intervals that are plotted on each frame, considering the set number of trail-values.
            
            % TODO (what? 2017-03-05)
            this.p_n_trail = fix(this.p_n_time_steps * this.p_trail_length);
            this.d_t_trails = cell(1,this.p_n_time_steps);
            for ii = 1:this.p_n_time_steps
                current_indexes = max(1,(ii-this.p_n_trail)*this.p_oversampling_factor) : 1 : ii*this.p_oversampling_factor;
                this.d_t_trails{ii} = current_indexes;
            end
        end
        
        % ---------------
        % DRAW-Methods
        % ---------------
        
        function [] = draw_init_gui_statics(this)
            % Prepares the positions of the figure and its two subplots.
            
            % TODO (what? 2017-03-05)
            
            this.h_fig.Position =   this.w_fig_position;
            
            this.h_sub1.Position =  this.w_sub1_position;
            this.h_sub2.Position =  this.w_sub2_position;
        end
        
        function [] = draw_init_gui_text_objects(this)
            % Prepares the annotations below the plots and populates the corresponding handle arrays.
            
            % Delete all previous instances of text boxes first as there
            % might be a different number of poles and zeros this time.
            for ii = 1:length(this.h_text_z_annot)
                this.h_text_z_annot(ii).delete
            end
            for ii = 1:length(this.h_text_p_annot)
                this.h_text_p_annot(ii).delete
            end
            
            % Text boxes for zero and pole contributions.
            annotation(this.h_fig,'TextBox',[this.w_border_horizontal_frac, (this.w_annotation_start_frac - 0.5*this.w_annotation_textbox_frac), this.w_plot_width_frac, this.w_annotation_textbox_frac],'String','Contribution of the zeros:','LineStyle','none','FontSize',9);
            for ii = 1:this.d_n_zeros
                this.h_text_z_annot(ii) = annotation(this.h_fig,'TextBox',[this.w_border_horizontal_frac, (this.w_annotation_start_frac - this.w_annotation_textbox_frac*(ii+1)), this.w_plot_width_frac, this.w_annotation_textbox_frac],'String',['zero ',num2str(ii)],'LineStyle','none','FontSize',9);
            end
            annotation(this.h_fig,'TextBox',[0.25, (this.w_annotation_start_frac - 0.5*this.w_annotation_textbox_frac), this.w_plot_width_frac, this.w_annotation_textbox_frac],'String','Contribution of the poles:','LineStyle','none','FontSize',9);
            for ii = 1:this.d_n_poles
                this.h_text_p_annot(ii) = annotation(this.h_fig,'TextBox',[0.25, (this.w_annotation_start_frac - this.w_annotation_textbox_frac*(ii+1)), this.w_plot_width_frac, this.w_annotation_textbox_frac],'String',['pole ',num2str(ii)],'LineStyle','none','FontSize',9);
            end

            % Text boxes for the results.
            annotation(this.h_fig,'TextBox',[0.5 max(0,(this.w_annotation_start_frac - 0.5*this.w_annotation_textbox_frac)) this.w_plot_width_frac this.w_annotation_textbox_frac],'String','Resulting value of G:','LineStyle','none','FontSize',9);
            this.h_text_res_annot(1) = annotation(this.h_fig,'TextBox',[0.5 max(0,(this.w_annotation_start_frac - 2*this.w_annotation_textbox_frac)) this.w_plot_width_frac this.w_annotation_textbox_frac],'String','resultline 1','LineStyle','none','FontSize',9);
            this.h_text_res_annot(2) = annotation(this.h_fig,'TextBox',[0.5 max(0,(this.w_annotation_start_frac - 3*this.w_annotation_textbox_frac)) this.w_plot_width_frac this.w_annotation_textbox_frac],'String','resultline 2','LineStyle','none','FontSize',9);
            
        end
        
        function [] = draw_init_plot_axes(this)
            % Prepares the geometry and behavior of the two subplots and attaches the proper callback methods to them.
            
            % TODO maybe separate into two separate functions for z and w
            % plot
            
            % z-plot ------------
            
            % Establish 1:1 aspect ratio.
            axis(this.h_sub1, 'equal');
            
            % Do this to let the user zoom all the way out with a double
            % click.
            xlim(this.h_sub1, 'manual'), ylim(this.h_sub1, 'manual');
            xlim(this.h_sub1, anp_stretch_centered(this.p_z_xspan,1.05)), ylim(this.h_sub1, anp_stretch_centered(this.p_z_yspan,1.05));
            zoom(this.h_fig, 'reset');
            xlim(this.h_sub1, this.p_z_xlim), ylim(this.h_sub1, this.p_z_ylim);
            
            % Default behavior of Matlab is to put the axes on the left and
            % bottom edge of the graph. We don't want this. Instead, place
            % the axes such that they go through the origin.
            this.h_sub1.XAxisLocation = 'origin';
            this.h_sub1.YAxisLocation = 'origin';

            grid(this.h_sub1, 'on');
            % -------------------
            
            % w-plot ------------
            % Same story as in z-plot.
            axis(this.h_sub2, 'equal');
            xlim(this.h_sub2, 'manual'), ylim(this.h_sub2, 'manual');
            xlim(this.h_sub2, anp_stretch_centered(this.p_w_xspan,1.05)), ylim(this.h_sub2, anp_stretch_centered(this.p_w_yspan,1.05));
            zoom(this.h_fig, 'reset');
            xlim(this.h_sub2, this.p_w_xlim), ylim(this.h_sub2, this.p_w_ylim);
            this.h_sub2.XAxisLocation = 'origin';
            this.h_sub2.YAxisLocation = 'origin';
            grid(this.h_sub2, 'on');
            % -------------------
            
            % Create zoom and pan objects and give them the handles to our
            % custom callback method.
            this.h_zoom = zoom(this.h_fig);
            this.h_zoom.ActionPostCallback = @this.cb_after_zoom_or_pan;
            this.h_pan = pan(this.h_fig);
            this.h_pan.ActionPostCallback = @this.cb_after_zoom_or_pan;
        end
        
        function [] = draw_init_line_plots(this)
            % Instantiates the full- and trail plot objects and sets their properties.
            
            % Prepare the full input- and output curves with some
            % transparency
            this.h_z_plot_full =            plot(this.h_sub1,0,0);
            set(this.h_z_plot_full,'Color',[0.05 0.4970 0.7410]);
            this.h_z_plot_full.Color(4) =   0.3;
            
            this.h_w_plot_full =            plot(this.h_sub2,0,0);
            set(this.h_w_plot_full,'Color',[0.05 0.4970 0.7410]);
            this.h_w_plot_full.Color(4) =   0.3;

            % Prepare the curves' trail plots and remember their handle for
            % later use.
            this.h_z_plot_trail =           plot(this.h_sub1,0,0);
            set(this.h_z_plot_trail,'Color',[255 215 0]/255,'linewidth',2);
            this.h_w_plot_trail =           plot(this.h_sub2,0,0);
            set(this.h_w_plot_trail,'Color',[255 215 0]/255,'linewidth',2);
        end
        
        function [] = draw_init_plot_arrows(this)
            % Prepare the arrow annotation objects used for the tip of the trails.
            
            % We're talking about the arrow annotations that mark the value
            % of the head of the trail at the current frame and remember
            % their handle for later use.
            this.h_z_plot_arrow =   annotation(this.h_fig, 'Arrow',[0,0.1],[0,0.1]);
            this.h_w_plot_arrow =   annotation(this.h_fig, 'Arrow',[1,0.9],[0,0.1]);
        end
        
        function [] = draw_init_z_plot_poles_zeros(this)
            % Draws x's at the location of poles and o's at the location of zeros in the left plot.
            
            % TODO we could be way more efficient if we knew which objects
            % we need to touch and which we don't. For the time being just
            % redo everything every time this method runs.
                        
            % First delete all existing objects to start fresh.
            for ii = 1:length(this.h_z_pz_objcts)
                this.h_z_pz_objcts{ii}.delete();
            end
            
            % Preallocate the cell-array containing the object handles.
            this.h_z_pz_objcts = cell(1,this.d_n_poles + this.d_n_zeros);
                        
            % Where's the center of the plot with the current limits?
            plot_center = mean(xlim(this.h_sub1)) + 1i*mean(ylim(this.h_sub1));
            
            objct_ii = 1;
            for p_ii = 1:this.d_n_poles
                current_pole = this.d_poles(p_ii);
                pole_trunc = this.trunc(current_pole, this.p_z_xlim, this.p_z_ylim);
                if current_pole ~= pole_trunc
                    % The current pole being not equal to its truncated value
                    % means that it lies outside the current subplot limits.
                    % Instead of drawing just an X at the border, plot an
                    % arrow, indicating that there's an outlier.
                    this.h_z_pz_objcts{objct_ii} = this.draw_text_arrow(this.h_sub1,[real(pole_trunc),imag(pole_trunc)],angle(current_pole - plot_center),this.p_z_arrow_length,' x',[255 140 0]/255);
                else
                    % The current pole lies within xlim and ylim, so draw
                    % an X at its location.
                    this.h_z_pz_objcts{objct_ii} = scatter(this.h_sub1,real(current_pole),imag(current_pole),60,'x','MarkerEdgeColor',[255 140 0]/255,'LineWidth',1.5);
                end
                objct_ii = objct_ii + 1;
            end
            
            % Same procedure as with poles but instead of X's, plot circles
            for z_ii = 1:this.d_n_zeros
                current_zero = this.d_zeros(z_ii);
                zero_trunc = this.trunc(current_zero, this.p_z_xlim, this.p_z_ylim);
                if current_zero ~= zero_trunc
                    this.h_z_pz_objcts{objct_ii} = this.draw_text_arrow(this.h_sub1,[real(zero_trunc),imag(zero_trunc)],angle(current_zero - plot_center),this.p_z_arrow_length,' o',[95 158 160]/255);
                else
                    this.h_z_pz_objcts{objct_ii} = scatter(this.h_sub1,real(current_zero),imag(current_zero),60,'o','MarkerEdgeColor',[70 130 180]/255,'LineWidth',1.5);
                end
                objct_ii = objct_ii + 1;
            end
        end
        
        % -----------------------------------------------------------------
        % GUI methods
        % -----------------------------------------------------------------
        function [] = ui_control_enable(this)
            % Makes all the custom buttons clickable.
            
            this.ui_btn_prev.Enable = 'on';
            this.ui_sw_pause.Enable = 'on';
            this.ui_btn_next.Enable = 'on';
            
            for ii = 1:length(this.ui_run_switches)
                this.ui_run_switches(ii).Enable = 'on';
            end

        end
        
        % -----------------------------------------------------------------
        % Methods for the running plot
        % -----------------------------------------------------------------
        
        function [] = draw_run_continuous_animation(this)
            % Manages the periodic update of the trails and their arrows.
            % Gets called whenever the animation should be turned on.
            
            % Whether an instance of this method already runs should be
            % checked prior to calling the method by testing for
            % this.s_draw_busy!
            
            try
                figure(this.h_fig);
                
                % Prevents data to be changed and is an indicator whether
                % this method already runs.
                this.s_draw_busy = true;
                
                tools.dbg('anp_gui[draw_run_continuous_animation]:\tStarting animation\n');
                
                % Main GUI update loop.
                while this.s_draw_allowed
                    
                    % Check whether the plot limits have been changed in
                    % the meantime during the last update below.
                    if this.s_check_limits
                        this.draw_update_limits_and_plots();
                        this.calc_z_arrow_length();
                        this.calc_w_arrow_length();
                        this.draw_init_z_plot_poles_zeros();
                        this.draw_one_frame();
                        this.s_check_limits = false;
                    end
                    
                    % Set the animation iterator to the next value.
                    this.a_time_ii = tools.iterator_modulo(this.a_time_ii + this.a_direction,this.p_n_time_steps);
                    
                    this.draw_one_frame();
                    
                    tools.dbg('anp_gui[draw_run_continuous_animation]:\tdrawing %d :-)\n',this.a_time_ii);
                    pause(1/10);
                end
                tools.dbg('anp_gui[draw_run_continuous_animation]:\tStopping animation\n');
                
                this.s_draw_busy = false;
            catch err
                % We ignore the exception when the figure is closed
                % in mid-update. Rethrow everything else.
                if ~strcmp(err.message,'Invalid or deleted object.')
                    rethrow(err);
                end
            end
        end
        
        function [] = draw_one_frame(this)
            % Manages the drawing tasks that one time step needs.
            
            this.draw_update_trails();
            this.draw_update_trail_head_arrows();
            this.draw_update_plot_titles();
            this.draw_update_textboxes();
        end
        
        function [] = draw_update_trails(this)
            % Updates the trail plots.
            
            set(this.h_z_plot_trail,'XData',real(this.d_z_values_truncated(this.d_t_trails{this.a_time_ii})),'YData',imag(this.d_z_values_truncated(this.d_t_trails{this.a_time_ii})));
            set(this.h_w_plot_trail,'XData',real(this.d_w_values_truncated(this.d_t_trails{this.a_time_ii})),'YData',imag(this.d_w_values_truncated(this.d_t_trails{this.a_time_ii})));        
        end
        
        function [] = draw_update_trail_head_arrows(this)
            % Draws the (arrow-)heads of the trails.
                        
            % The trail arrows need to point always into the same
            % direction, no matter in which direction the animation goes.
            % So we get the "previous" and current head location and draw
            % the arrow at the current location with the direction from the
            % previous one.
            current_values_index =  this.a_time_ii * this.p_oversampling_factor;
            prev_values_index =     tools.iterator_modulo(this.a_time_ii * this.p_oversampling_factor - 1,this.p_n_time_steps * this.p_oversampling_factor);
            z_values_head_prev =    this.d_z_values(prev_values_index);
            w_values_head_prev =    this.d_w_values(prev_values_index);
            z_values_head =         this.d_z_values(current_values_index);
            w_values_head =         this.d_w_values(current_values_index);
            
            
            z_phi =                 angle(z_values_head - z_values_head_prev);
            w_phi =                 angle(w_values_head - w_values_head_prev);
            
            
            z_arrow_tip_x = real(this.d_z_values_truncated(this.a_time_ii * this.p_oversampling_factor));
            z_arrow_tip_y = imag(this.d_z_values_truncated(this.a_time_ii * this.p_oversampling_factor));
            
            this.draw_update_arrow(this.h_sub1,this.h_z_plot_arrow,z_arrow_tip_x,z_arrow_tip_y,z_phi,this.p_z_arrow_length);
            
            
            w_arrow_tip_x = real(this.d_w_values_truncated(this.a_time_ii * this.p_oversampling_factor));
            w_arrow_tip_y = imag(this.d_w_values_truncated(this.a_time_ii * this.p_oversampling_factor));
            
            this.draw_update_arrow(this.h_sub2,this.h_w_plot_arrow,w_arrow_tip_x,w_arrow_tip_y,w_phi,this.p_w_arrow_length);
        end
        
        function [] = draw_update_plot_titles(this)
            % Updates the text of the plot titles.
            % We calculate the angles directly from the value at the arrow
            % tip, so the phase can be off by multiples of 360� compared to
            % the cumulative phase in the text boxes below the graph.
            
            current_z_value = this.d_z_values(this.a_time_ii * this.p_oversampling_factor);
            current_w_value = this.d_w_values(this.a_time_ii * this.p_oversampling_factor);
            
            title(this.h_sub1,['Current value of D-Curve: ',num2str(current_z_value,'%.1f'),': M = ',num2str(abs(current_z_value),'%.2f'),' p = ',num2str(rad2deg(angle(current_z_value)),'%.2f'),'�']);
            title(this.h_sub2,['Nyquist: G(',num2str(current_z_value,'%.1f'),') = ',num2str(current_w_value,'%.1f'),': M = ',num2str(abs(current_w_value),'%.2f'),' p = ',num2str(rad2deg(angle(current_w_value)),'%.2f'),'�']);
        end
        
        function [] = draw_update_textboxes(this)
            % Updates the zero- and pole contribution and the cumulative values and puts the information in the pre-allocated textboxes.
                        
            res_magnitude = 1;
            res_phase =     0;
            
            % *************************
            % Text before contributions
            % *************************
            
            switch this.d_n_zeros
                case 0
                    res_magnitude_txt =         '1';
                case 1
                    res_magnitude_txt =         '';
                otherwise
                    res_magnitude_txt =         '[';
            end
            res_phase_txt =                     '[';
            
            % ******************
            % Zero contributions
            % ******************
            
            for z = 1:this.d_n_zeros
                % Calculate (s-z_i).
                z_contribution =                this.d_z_values(this.a_time_ii * this.p_oversampling_factor) - this.d_zeros(z);
                res_magnitude =                 res_magnitude * abs(z_contribution);
                res_phase =                     res_phase + angle(z_contribution);
                
                % Put together the text for the text box that describes the
                % current zero's contribution, for example
                % Z1: (0.0+1.1i) - (-0.7): M=1.3 p=57.9�
                this.h_text_z_annot(z).String = ['Z',num2str(z),    ': (',  num2str(this.d_z_values(this.a_time_ii * this.p_oversampling_factor),           '%.1f'),') - (',num2str(this.d_zeros(z),'%.1f'),'): M=',num2str(abs(z_contribution),'%.1f'),' p=',num2str(rad2deg(angle(z_contribution)),'%.1f'),'�'];
                if z == 1
                    res_magnitude_txt =      	[res_magnitude_txt,         num2str(abs(z_contribution),           '%.2f')     ];
                else
                    res_magnitude_txt =        	[res_magnitude_txt, ' * ',  num2str(abs(z_contribution),           '%.2f')     ];
                end
                res_phase_txt =             [res_phase_txt,     ' + (', num2str(rad2deg(angle(z_contribution)),'%.2f'), ')'];
            end
            
            % ****************************************
            % Text between zero and pole contributions
            % ****************************************
            
            if this.d_n_zeros <= 1
                res_magnitude_txt =           	[res_magnitude_txt, ' / '];
            else
                res_magnitude_txt =            	[res_magnitude_txt, '] / '];
            end
            
            if this.d_n_poles >= 2
                res_magnitude_txt =            	[res_magnitude_txt, '['];
            end
            
            res_phase_txt =                     [res_phase_txt,     '] ['];

            
            % ******************
            % Pole contributions
            % ******************
            
            % Same story as with the zeros.
            for p = 1:this.d_n_poles
                p_contribution =                this.d_z_values(this.a_time_ii * this.p_oversampling_factor) - this.d_poles(p);
                res_magnitude =                 res_magnitude / abs(p_contribution);
                res_phase =                     res_phase - angle(p_contribution);
                
                this.h_text_p_annot(p).String = ['P',num2str(p),    ': (',  num2str(this.d_z_values(this.a_time_ii * this.p_oversampling_factor),           '%.1f'),') - (',num2str(this.d_poles(p),'%.1f'),'): M=',num2str(abs(p_contribution),'%.1f'),' p=',num2str(rad2deg(angle(p_contribution)),'%.1f'),'�'];            
                if p == 1
                    res_magnitude_txt =     	[res_magnitude_txt,         num2str(abs(p_contribution),           '%.2f')     ];
                else
                    res_magnitude_txt =      	[res_magnitude_txt, ' * ',  num2str(abs(p_contribution),           '%.2f')     ];
                end
                res_phase_txt =            	[res_phase_txt,     ' - (', num2str(rad2deg(angle(p_contribution)),'%.2f'), ')'];
            end
            
            % ************************
            % Text after contributions
            % ************************
            
            if this.d_n_poles >= 2
                res_magnitude_txt =             [res_magnitude_txt,']'];
            end
            res_phase_txt =                   	[res_phase_txt,    ']'];
            
            res_magnitude_txt =               	[res_magnitude_txt,' = ',   num2str(res_magnitude,      '%.3f')];
            res_phase_txt =                    	[res_phase_txt,    ' = ',   num2str(rad2deg(res_phase),'%.3f'),'�'];
            
            this.h_text_res_annot(1).String =   ['Magnitude: ',  res_magnitude_txt];
            this.h_text_res_annot(2).String =   ['Phase:       ',res_phase_txt];
        end
        
        function [] = draw_update_limits_and_plots(this)
            % Updates the axes' limits and their plot data.
            % Gets called by GUI callback-methods. For example after a zoom
            % or pan event.
            
            this.p_z_xlim = xlim(this.h_sub1);
            this.p_z_ylim = ylim(this.h_sub1);

            this.calc_truncated_z_values();
            this.draw_update_full_z_plot();

            this.p_w_xlim = xlim(this.h_sub2);
            this.p_w_ylim = ylim(this.h_sub2);
            
            this.calc_truncated_w_values();
            this.draw_update_full_w_plot();
            
            tools.dbg('anp_gui[draw_update_limits_and_plots]:\tz_xlim=[%.2f,%.2f], z_ylim =[%.2f,%.2f], w_xlim=[%.2f,%.2f], w_ylim=[%.2f,%.2f]\n',this.p_z_xlim(1),this.p_z_xlim(2),this.p_z_ylim(1),this.p_z_ylim(2),this.p_w_xlim(1),this.p_w_xlim(2),this.p_w_ylim(1),this.p_w_ylim(2));
        end
        
        function [] = draw_update_full_z_plot(this)
            % Wrapper method. Updates the z-plot's data.
            
            set(this.h_z_plot_full,'XData',real(this.d_z_values_truncated),'YData',imag(this.d_z_values_truncated));
        end
        
        function [] = draw_update_full_w_plot(this)
            % Wrapper method. Updates the w-plot's data.
            set(this.h_w_plot_full,'XData',real(this.d_w_values_truncated),'YData',imag(this.d_w_values_truncated));
        end
        
        
        % -----------------------------------------------------------------
        % Callback methods
        % -----------------------------------------------------------------
        
        function [] = cb_step(this,~,~,dir) % ignored parameters are src and evt
            % CB method reacting to a "step DIR" button.
            % The dir value is hardcoded in the button's instantiation.
            
            this.s_draw_busy = true;
            tools.dbg('anp_gui[cb_step]:\tStepping button.\n');
            
            this.a_time_ii = tools.iterator_modulo(this.a_time_ii + dir,this.p_n_time_steps);
            this.draw_one_frame();
            
            this.s_draw_busy = false;
        end
        
        
        function [] = cb_run(this,src,~,dir) % ignored parameter is evt
            % CB method reacting to one of the "play" buttons.
            % The dir values are hardcoded in the buttons' instantiations.
            
            if this.s_draw_busy && (dir == this.a_direction)
                % If draw_run_continuous_animation() already runs with the
                % same dir as the requested one, do nothing.
                
                return;
            elseif this.s_draw_busy && (dir ~= this.a_direction)
                % If draw_run_continuous_animation() already runs, but with
                % a different dir, change it.
                
                this.a_direction = dir;
                
                % This case must be because of another "play" button being
                % pressed, so turn off the old one.
                for ii = 1:length(this.ui_run_switches)
                    if ~isequal(src,this.ui_run_switches(ii))
                        this.ui_run_switches(ii).State = 'off';
                    end
                end
            else
                % draw_run_continuous_animation() wasn't running yet.
                
                tools.dbg('anp_gui[cb_run]:\tStart requested.\n');
                if this.s_gui_nominal
                    % The GUI and data are ready
                    
                    this.ui_sw_pause.State =    'off';
                    
                    % Block the "step" buttons while the animation
                    % is running.
                    this.ui_btn_prev.Enable =   'off';
                    this.ui_btn_next.Enable =   'off';
                    
                    % Turn off the other "play" buttons.
                    for ii = 1:length(this.ui_run_switches)
                        if ~isequal(src,this.ui_run_switches(ii))
                            this.ui_run_switches(ii).State = 'off';
                        end
                    end
                    
                    % Set the animation direction to the requested one.
                    this.a_direction =          dir;
                    
                    % s_draw_allowed is a means to stop
                    % draw_run_continuous_animation(), so better set
                    % it to true.
                    this.s_draw_allowed =       true;
                    
                    this.draw_run_continuous_animation();
                else
                    warning('Data not ready!');
                    src.State =  'off';
                end
            end
        end
        
        function [] = cb_run_switch_off(this,src,~,dir)
            % CB that doesn't let a "play" button to be turned off.
            if this.s_draw_allowed == true && (dir == this.a_direction)
                src.State = 'on';
            end
        end
        
        function [] = cb_pause(this,src,~)
            % CB method that reacts to the "pause" button.
            
            tools.dbg('anp_gui[cb_pause]:\tPause requested.\n');
            
            % This causes draw_run_continuous_animation() to stop.
            this.s_draw_allowed =   false;
            
            % Set the "play" and "pause" button states correctly.
            src.State =             'on';
            for ii = 1:length(this.ui_run_switches)
                this.ui_run_switches(ii).State = 'off';
            end
            
            % Unblock the "step" buttons.
            this.ui_btn_prev.Enable =   'on';
            this.ui_btn_next.Enable =   'on';
        end
        
        function cb_after_zoom_or_pan(this,~,~)
            % CB methods for zoom and pan events.
            
            if ~this.s_draw_busy
                % As long as draw_run_continuous_animation() isn't running,
                % we do the updates directly.
                
                this.draw_update_limits_and_plots();
                this.calc_z_arrow_length();
                this.calc_w_arrow_length();
                this.draw_init_z_plot_poles_zeros();
                this.draw_one_frame();
            else
                % draw_run_continuous_animation() is still running, notify
                % it of the changed axes limits.
                
                this.s_check_limits =   true;
            end
        end
        
        % ---------------------
        % Utility functions
        % ---------------------
        
        % UNUSED
        function arrowHandle = draw_arrow(this,parent,x0,phi,l)
            % Places an arrow annotation at the given location.
            r = l*[cos(phi),sin(phi)];
            
            arrowHandle = annotation(this.h_fig,'Arrow');
            set(arrowHandle,'parent',parent,'position',[x0-r,r]);
        end
        
        function h_arrow = draw_text_arrow(this,parent,x0,phi,l,text,color)
            % Places a text arrow annotation at the given location.
            r = l*[cos(phi),sin(phi)];

            h_arrow = annotation(this.h_fig,'TextArrow');
            set(h_arrow,'parent',parent,'position',[x0-r,r],'String',text,'Color',color);
        end
        
        function [] = draw_update_arrow(~,h_plot,h_arrow,x,y,phi,length) % ignored parameter is 'this'
            % Updates an annotation('arrow')-object's position
            
            % First calculate the start and end coordinates within the plot
            dx = length * cos(phi);
            dy = length * sin(phi);
            
            arrow_end_x = x - dx;
            arrow_end_y = y - dy;
            
            % Transform the plot coordinates to relative position
            % coordinates in the figure, as this is the arrow's parent
            % object.
            [transformed_x,transformed_y] = ds2nfu(h_plot,[arrow_end_x,x],[arrow_end_y,y]);
            
            h_arrow.X = transformed_x;
            h_arrow.Y = transformed_y;
        end

        function values_truncated = trunc(~,values,xlim,ylim) % ignored parameter is 'this'
            % Clips the input values to the maximum value allowed inside the plot.
            
            % We actually don't clip to the maximum value but only 97% so
            % the clipped value is well inside the plot.
            xlim = anp_stretch_centered(xlim,0.97);
            ylim = anp_stretch_centered(ylim,0.97);
            values_truncated = max(xlim(1), min(xlim(2), real(values))) ...
                          + 1i*max(ylim(1), min(ylim(2), imag(values)));
        end
    end
end