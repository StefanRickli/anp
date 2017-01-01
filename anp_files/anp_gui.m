% as in p.41
% last line: 516
classdef anp_gui < handle
    properties(SetAccess = private)
        % g: general properties
        g_uid               % Java String
        g_legacy            string
        
        % s: state
        s_data_ready        logical
        s_draw_allowed      logical
        s_draw_busy         logical
        s_check_limits      logical
        
        % h: handles
        h_anp_tf_processor
        
        h_fig               matlab.ui.Figure
        h_sub1              matlab.graphics.axis.Axes
        h_sub2              matlab.graphics.axis.Axes
        h_zoom              matlab.graphics.interaction.internal.zoom
        h_pan               matlab.graphics.interaction.internal.pan
        
        h_text_z_annot      matlab.graphics.shape.TextBox
        h_text_p_annot      matlab.graphics.shape.TextBox
        h_text_res_annot    matlab.graphics.shape.TextBox
        
        h_z_plot_full       matlab.graphics.chart.primitive.Line
        h_z_plot_trail      matlab.graphics.chart.primitive.Line
        h_z_plot_arrow      matlab.graphics.shape.Arrow
        h_z_pz_objcts       % array of handles with possibly mixed type
        
        h_w_plot_full       matlab.graphics.chart.primitive.Line
        h_w_plot_trail      matlab.graphics.chart.primitive.Line
        h_w_plot_arrow      matlab.graphics.shape.Arrow
                
        % ui: GUI elements
        ui_icons            double
        ui_toolbar          matlab.ui.container.Toolbar
        ui_run_switches     matlab.ui.container.toolbar.ToggleTool
        ui_sw_pause         matlab.ui.container.toolbar.ToggleTool
        ui_btn_prev         matlab.ui.container.toolbar.PushTool
        ui_btn_next         matlab.ui.container.toolbar.PushTool
        
        % a: animation
        a_direction         double
        a_time_ii           double

        % w: window properties
        w_fig_position      double
        w_sub1_position     double
        w_sub2_position     double
        
        w_border                    double
        w_plot_size                 double
        w_plot_width_frac           double
        w_border_horizontal_frac    double
        w_border_vertical_frac      double
        w_annotation_start_frac     double
        w_annotation_textbox_frac   double
        
        % p: plot properties (z: left, w=f(z): right)
        p_n_time_steps      double
        p_oversampling_factor double
        n_data_points       double
        p_trail_length      double
        p_n_trail           double
        
        p_z_x0              double
        p_z_dims            double
        p_z_auto_lims       logical
        p_z_xlim            double
        p_z_xspan           double
        p_z_ylim            double
        p_z_yspan           double
        p_z_width           double
        p_z_height          double
        p_z_arrow_length    double
        
        p_w_x0              double
        p_w_dims            double
        p_w_auto_lims       logical
        p_w_xlim            double
        p_w_xspan           double
        p_w_ylim            double
        p_w_yspan           double
        p_w_width           double
        p_w_height          double
        p_w_arrow_length    double
        
        % d: data
        d_n_poles               double
        d_n_zeros               double
        d_poles                 double
        d_zeros                 double
        d_R                     double
        
        d_t_values              double
        d_t_oversampled         double
        d_t_trails
        
        d_z_values              double
        d_z_values_truncated    double
        
        d_w_values              double
        d_w_values_truncated    double        
    end
    methods(Access = public)
        function this = anp_gui()
            this.g_uid =        java.util.UUID.randomUUID.toString;
            this.g_legacy =     anp_check_Matlab_version();
            this.s_data_ready = false;
            this.s_draw_allowed = false;
            this.s_draw_busy =  false;
            this.a_time_ii =    1;
            this.a_direction =  1;
            
            this.h_fig =        figure('CloseRequestFcn',@this.on_figure_delete);
            this.h_sub1 =       subplot(1,2,1);
            hold on;
            this.h_sub2 =       subplot(1,2,2);
            hold on;
            
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
            % TODO check for dirty values and later redraw
            
            this.set_z_plot_limits(new_props);
            this.set_w_plot_limits(new_props);
            this.set_poles_zeros(new_props);
            this.set_window_props(new_props);
            this.set_plot_props(new_props);
            
            this.set_anp_tf_processor(new_props);
        end
        
        function [] = set_window_props(this,new_props)
            this.w_border =     new_props.border;
            this.w_plot_size =  new_props.plot_size;
        end
        
        function [] = set_plot_props(this,new_props)
            this.p_trail_length =     new_props.trail_length;
        end
        
        function [] = set_z_plot_limits(this,new_props)
            % TODO check for proper array dimensions
            
            this.p_z_x0 =           new_props.z_plot_x0;
            this.p_z_dims =         new_props.z_plot_dims;
            this.p_z_auto_lims =    new_props.z_plot_auto_lims;
        end
        
        function [] = set_w_plot_limits(this,new_props)
            % TODO check for proper array dimensions
            
            this.p_w_x0 =           new_props.w_plot_x0;
            this.p_w_dims =         new_props.w_plot_dims;
            this.p_w_auto_lims =    new_props.w_plot_auto_lims;
        end
        
        function [] = set_radii(this,new_props)
            % TODO check for proper array dimensions
            % TODO radii comes from tf_processor!
            this.d_radii =          new_props.radii;
        end
                
        function [] = set_poles_zeros(this,new_props)
            % TODO check for proper array dimensions
            
            this.d_poles =          new_props.tf_poles;
            this.d_zeros =          new_props.tf_zeros;
            
            this.d_n_poles =        length(new_props.tf_poles);
            this.d_n_zeros =        length(new_props.tf_zeros);
        end

        function [] = set_anp_tf_processor(this,new_props)
            this.h_anp_tf_processor = new_props.processor_handle;
        end
        
        function [] = init_visuals(this)
            this.fetch_R();
            this.fetch_data();
            this.calc_gui_positions();
            this.calc_plot_axis_limits();
            this.calc_plot_z_arrow_length();
            this.calc_plot_w_arrow_length();
            this.calc_truncated_z_values();
            this.calc_truncated_w_values();
            this.calc_trail_indexes();
            this.draw_init_gui_statics();
            this.draw_init_gui_text_objects();
            this.draw_init_plot_axes();
            this.draw_init_line_plots();
            this.draw_init_plot_arrows();
            this.draw_init_z_plot_poles_zeros();
            this.draw_update_full_z_plot();
            this.draw_update_full_w_plot();
            this.draw_one_frame();
            
            this.s_data_ready =     true;
            this.ui_control_enable();
        end
        
        
        function delete(this)
            tools.dbg('anp_gui[delete]\t%s: Deletion requested.\n',this.g_uid);
            
            if isvalid(this.h_fig)
                close(this.h_fig);
            end
        end
    end
    methods(Access = private)
        function on_figure_delete(this,src,~) % last argument is 'evt'
            % TODO remove all event listeners here!
            tools.dbg('anp_gui[on_figure_delete]:\tDeleting figure window.\n');
            
            delete(src);
            this.delete();
        end
        
        function [] = load_icons(this)
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
            assert(isvalid(this.h_anp_tf_processor));
            
            this.d_R =                  this.h_anp_tf_processor.get_R();
        end
        
        function [] = fetch_data(this)
            assert(isvalid(this.h_anp_tf_processor));
            
            time_data =                     this.h_anp_tf_processor.get_time_points();
            this.d_t_values =               time_data.time_points;
            this.d_t_oversampled =          time_data.data_points;
            this.p_n_time_steps =           time_data.time_props.n_time_steps;
            this.p_oversampling_factor =    time_data.time_props.oversampling_factor;
            this.n_data_points =            time_data.time_props.n_data_points;
            
            function_data = this.h_anp_tf_processor.get_data();
            this.d_z_values =       function_data.z_values;
            this.d_w_values =       function_data.w_values;
        end
        
        function [] = calc_gui_positions(this)
            % calculate the correct positions (relative to interior of the figure)
            fig_plot_height = 2*this.w_border + this.w_plot_size;              % plot + border above and below, pixel
            fig_annotation_textbox_height = 14;                                                 % one box, pixel
            fig_annotation_height_sum = (max(this.d_n_zeros,this.d_n_poles)+2)*fig_annotation_textbox_height; % cumulative, pixel
            switch(this.g_legacy)
                case 'R2015b_or_newer'
                    fig_legacy_correction = 0;
                otherwise
                    fig_legacy_correction = 3*fig_annotation_textbox_height;
            end
            fig_height = fig_plot_height + fig_annotation_height_sum + fig_legacy_correction;   % pixel
            fig_width = 3*this.w_border + 2*this.w_plot_size;                                  % pixel
            
            % fig.Position expects pixels as unit
            this.w_fig_position = [100 100 fig_width fig_height];
            
            % as subplot.Postion expects fractions of the inside of the figure, we
            % recalculate the pixel values
            this.w_border_horizontal_frac = this.w_border/fig_width;                      % fraction
            this.w_border_vertical_frac = this.w_border/fig_height;                       % fraction
            this.w_plot_width_frac = (fig_width - 3*this.w_border)/(2*fig_width);         % fraction
            fig_plot_height_frac = 1-(2*this.w_border + fig_annotation_height_sum)/fig_height; % fraction
            
            % subX.Position expects fractions of the inside of the figure as unit
            this.w_sub1_position = [this.w_border_horizontal_frac, (1 - this.w_border_vertical_frac - fig_plot_height_frac), this.w_plot_width_frac, fig_plot_height_frac];
            this.w_sub2_position = [(2*this.w_border_horizontal_frac + this.w_plot_width_frac),(1 - this.w_border_vertical_frac - fig_plot_height_frac), this.w_plot_width_frac, fig_plot_height_frac];
            
            % text box object parameters:
            % again we need to know relative (fractions) positions inside the figure
            this.w_annotation_start_frac = (fig_annotation_height_sum + 10)/fig_height;    % fraction
            this.w_annotation_textbox_frac = fig_annotation_textbox_height/fig_height;     % fraction
        end
        
        function [] = calc_plot_axis_limits(this)
            % TODO make sure that all props and data are set before calling
            % this method
            
            assert(~any([isempty([this.d_zeros,this.d_poles]),...
                         isempty(this.d_R) || isnan(this.d_R),...
                         isempty(this.d_w_values),...
                         isempty(this.p_z_auto_lims),...
                         isempty(this.p_w_auto_lims)]));
            
            if this.p_z_auto_lims
                % try to find optimal axis limits for the two plots, or let the user
                % decide if he chose to set either 'auto_size' to false
                
                [this.p_z_xlim,this.p_z_ylim] = anp_plot_auto_zoom_z([this.d_zeros,this.d_poles],this.d_R);
            else
               this.p_z_xlim = [(this.p_z_x0(1) - this.p_z_dims(1)/2),(this.p_z_x0(1) + this.p_z_dims(1)/2)];
               this.p_z_ylim = [(this.p_z_x0(2) - this.p_z_dims(2)/2),(this.p_z_x0(2) + this.p_z_dims(2)/2)];
            end
            
            [this.p_z_xspan,this.p_z_yspan] = anp_plot_find_span([this.d_zeros,this.d_poles,this.d_z_values]);
            if diff(this.p_z_xlim) < 100*eps
                this.p_z_xlim = this.p_z_xspan;
            end
            if diff(this.p_z_ylim) < 100*eps
                this.p_z_ylim = this.p_z_yspan;
            end
            
            
            if this.p_w_auto_lims
                [this.p_w_xlim,this.p_w_ylim] = anp_plot_auto_zoom_w(this.d_w_values);
            else
               this.p_w_xlim = [(this.p_w_x0(1) - this.p_w_dims(1)/2),(this.p_w_x0(1) + this.p_w_dims(1)/2)];
               this.p_w_ylim = [(this.p_w_x0(2) - this.p_w_dims(2)/2),(this.p_w_x0(2) + this.p_w_dims(2)/2)];
            end
            
            [this.p_w_xspan,this.p_w_yspan] = anp_plot_find_span(this.d_w_values);
            if diff(this.p_w_xlim) < 100*eps
                this.p_w_xlim = this.p_w_xspan;
            end
            if diff(this.p_w_xlim) < 100*eps
                this.p_w_ylim = this.p_w_yspan;
            end
        end
        
        function [] = calc_plot_z_arrow_length(this)
            in_axis_width = diff(this.p_z_xlim);        % [1]
            in_axis_height = diff(this.p_z_ylim);       % [1]
            this.p_z_arrow_length = 0.04 * sqrt(in_axis_width^2 + in_axis_height^2);      % [1]
            tools.dbg('anp_gui[calc_plot_z_arrow_length]:\t%.5f\n',this.p_z_arrow_length);
        end
        
        function [] = calc_plot_w_arrow_length(this)
            % TODO handle very short arrow lengths better
            % This involves setting the position relative to the figure/subplot
            % and not within the coordinate system, which is cumbersome...
            out_axis_width = diff(this.p_w_xlim);     % [1]
            out_axis_height = diff(this.p_w_ylim);    % [1]
            this.p_w_arrow_length = 0.04 * sqrt(out_axis_width^2 + out_axis_height^2);   % [1]
            if this.p_w_arrow_length < 0.05
                this.p_w_arrow_length = 0.0001;
            end
            tools.dbg('anp_gui[calc_plot_w_arrow_length]:\t%.5f\n',this.p_w_arrow_length);
        end
        
        function [] = calc_truncated_z_values(this)
            this.d_z_values_truncated =  this.trunc(this.d_z_values,this.p_z_xlim,this.p_z_ylim);
        end
        function [] = calc_truncated_w_values(this)
            this.d_w_values_truncated =  this.trunc(this.d_w_values,this.p_w_xlim,this.p_w_ylim);
        end
        
        % prepare t-intervals that are plotted on each frame, considering the
        % set number of trail-values
        function [] = calc_trail_indexes(this) % TODO
            this.p_n_trail = fix(this.p_n_time_steps * this.p_trail_length);
            this.d_t_trails = cell(1,this.p_n_time_steps);
            for ii = 1:this.p_n_time_steps
                this.d_t_trails{ii} = max(1,(ii-this.p_n_trail)*this.p_oversampling_factor):ii*this.p_oversampling_factor;
            end
        end
        
        function [] = draw_init_gui_statics(this)
            % TODO 
            this.h_fig.Position =   this.w_fig_position;
            
            this.h_sub1.Position =  this.w_sub1_position;
            this.h_sub2.Position =  this.w_sub2_position;
        end
        
        % prepare the annotations below the plots and remember their handles
        function [] = draw_init_gui_text_objects(this)
            % parameters get set in 'calc_gui_positions'
            
            for ii = 1:length(this.h_text_z_annot)
                this.h_text_z_annot(ii).delete
            end
            for ii = 1:length(this.h_text_p_annot)
                this.h_text_p_annot(ii).delete
            end
            
            % zero and pole contributions
            annotation(this.h_fig,'TextBox',[this.w_border_horizontal_frac, (this.w_annotation_start_frac - 0.5*this.w_annotation_textbox_frac), this.w_plot_width_frac, this.w_annotation_textbox_frac],'String','Contribution of the zeros:','LineStyle','none','FontSize',9);
            for ii = 1:this.d_n_zeros
                this.h_text_z_annot(ii) = annotation(this.h_fig,'TextBox',[this.w_border_horizontal_frac, (this.w_annotation_start_frac - this.w_annotation_textbox_frac*(ii+1)), this.w_plot_width_frac, this.w_annotation_textbox_frac],'String',['zero ',num2str(ii)],'LineStyle','none','FontSize',9);
            end
            annotation(this.h_fig,'TextBox',[0.25, (this.w_annotation_start_frac - 0.5*this.w_annotation_textbox_frac), this.w_plot_width_frac, this.w_annotation_textbox_frac],'String','Contribution of the poles:','LineStyle','none','FontSize',9);
            for ii = 1:this.d_n_poles
                this.h_text_p_annot(ii) = annotation(this.h_fig,'TextBox',[0.25, (this.w_annotation_start_frac - this.w_annotation_textbox_frac*(ii+1)), this.w_plot_width_frac, this.w_annotation_textbox_frac],'String',['pole ',num2str(ii)],'LineStyle','none','FontSize',9);
            end

            % cumulative calculations
            annotation(this.h_fig,'TextBox',[0.5 max(0,(this.w_annotation_start_frac - 0.5*this.w_annotation_textbox_frac)) this.w_plot_width_frac this.w_annotation_textbox_frac],'String','Resulting value of G:','LineStyle','none','FontSize',9);
            this.h_text_res_annot(1) = annotation(this.h_fig,'TextBox',[0.5 max(0,(this.w_annotation_start_frac - 2*this.w_annotation_textbox_frac)) this.w_plot_width_frac this.w_annotation_textbox_frac],'String','resultline 1','LineStyle','none','FontSize',9);
            this.h_text_res_annot(2) = annotation(this.h_fig,'TextBox',[0.5 max(0,(this.w_annotation_start_frac - 3*this.w_annotation_textbox_frac)) this.w_plot_width_frac this.w_annotation_textbox_frac],'String','resultline 2','LineStyle','none','FontSize',9);
            
        end
        
        function [] = draw_init_plot_axes(this)
            % TODO maybe separate into two separate functions for z and w
            % plot
            
            subplot(this.h_sub1);
            axis equal; % for 1:1 aspect ratio
            xlim manual, ylim manual;
            xlim(anp_stretch_centered(this.p_z_xspan,1.05)), ylim(anp_stretch_centered(this.p_z_yspan,1.05));
            zoom reset;
            xlim(this.p_z_xlim), ylim(this.p_z_ylim);
            anp_plot_axes_origin(this.g_legacy);
            grid on;
            
            subplot(this.h_sub2);
            axis equal;
            xlim manual, ylim manual;
            xlim(anp_stretch_centered(this.p_w_xspan,1.05)), ylim(anp_stretch_centered(this.p_w_yspan,1.05));
            zoom reset;
            xlim(this.p_w_xlim), ylim(this.p_w_ylim);
            anp_plot_axes_origin(this.g_legacy);
            grid on;
            
            this.h_zoom = zoom;
            this.h_zoom.ActionPostCallback = @this.cb_after_zoom_or_pan;
            this.h_pan = pan;
            this.h_pan.ActionPostCallback = @this.cb_after_zoom_or_pan;
        end
        
        function [] = draw_init_line_plots(this)
            % plot the full input- and output curves with some transparency
            this.h_z_plot_full =            plot(this.h_sub1,0,0);
            set(this.h_z_plot_full,'Color',[0.05 0.4970 0.7410]);
            this.h_z_plot_full.Color(4) =   0.3;
            
            this.h_w_plot_full =            plot(this.h_sub2,0,0);
            set(this.h_w_plot_full,'Color',[0.05 0.4970 0.7410]);
            this.h_w_plot_full.Color(4) =   0.3;

            % prepare the curves' trail plots and remember their handle for later
            % use
            this.h_z_plot_trail =           plot(this.h_sub1,0,0);
            set(this.h_z_plot_trail,'Color',[255 215 0]/255,'linewidth',2);
            this.h_w_plot_trail =           plot(this.h_sub2,0,0);
            set(this.h_w_plot_trail,'Color',[255 215 0]/255,'linewidth',2);
        end
        
        function [] = draw_init_plot_arrows(this)
            % prepare the arrows annotations that mark the value of the head of the
            % trail at the current frame and remember their handle for later use
            subplot(this.h_sub1);
            this.h_z_plot_arrow =       annotation('Arrow',[0 0],[1 0]);

            subplot(this.h_sub2);
            this.h_w_plot_arrow =       annotation('Arrow',[0 0],[1 0]);
        end
        
        % draw the this.d_poles and this.d_zeros in the left (input function) subplot
        function [] = draw_init_z_plot_poles_zeros(this)
            for ii = 1:length(this.h_z_pz_objcts)
                this.h_z_pz_objcts.delete();
            end
            this.h_z_pz_objcts = [];
            
            subplot(this.h_sub1);
            
            for p_ii = 1:this.d_n_poles
                current_pole = this.d_poles(p_ii);
                pole_trunc = this.trunc(current_pole, this.p_z_xlim, this.p_z_ylim);
                if current_pole ~= pole_trunc
                    % The current pole being not equal to its truncated value
                    % means that it lies outside the current subplot limits.
                    % Instead of drawing just a circle at the border, plot an
                    % arrow, indicating that there's an outlier.
                    this.h_z_pz_objcts = drawTextArrow([real(pole_trunc),imag(pole_trunc)],angle(current_pole),in_arrow_length,' x',[255 140 0]/255);
                else
                    this.h_z_pz_objcts = scatter(real(current_pole),imag(current_pole),60,'x','MarkerEdgeColor',[255 140 0]/255,'LineWidth',1.5);
                end
            end
            
            for z_ii = 1:this.d_n_zeros
                current_zero = this.d_zeros(z_ii);
                zero_trunc = this.trunc(current_zero, this.p_z_xlim, this.p_z_ylim);
                if current_zero ~= zero_trunc
                    this.h_z_pz_objcts = drawTextArrow([real(zero_trunc),imag(zero_trunc)],angle(current_zero),in_arrow_length,' o',[95 158 160]/255);
                else
                    this.h_z_pz_objcts = scatter(real(current_zero),imag(current_zero),60,'o','MarkerEdgeColor',[70 130 180]/255,'LineWidth',1.5);
                end
            end
        end
        
        % -----------------------------------------------------------------
        % GUI methods
        % -----------------------------------------------------------------
        function [] = ui_control_enable(this)
            this.ui_btn_prev.Enable = 'on';
            this.ui_sw_pause.Enable = 'on';
            this.ui_btn_next.Enable = 'on';
            
            for ii = 1:length(this.ui_run_switches)
                this.ui_run_switches(ii).Enable = 'on';
            end

        end
        
        % -----------------------------------------------------------------
        % methods for the running plot
        % -----------------------------------------------------------------
        function [] = draw_run_continuous_animation(this)
            try
                this.s_draw_busy = true;
                tools.dbg('anp_gui[draw_run_continuous_animation]:\tStarting animation\n');

                while this.s_draw_allowed
                    if this.s_check_limits
                        this.draw_update_limits_and_plots();
                        this.calc_plot_z_arrow_length();
                        this.calc_plot_w_arrow_length();
                        this.draw_one_frame();
                        this.s_check_limits = false;
                    end
                    
                    this.a_time_ii = tools.iterator_modulo(this.a_time_ii + this.a_direction,this.p_n_time_steps);

                    this.draw_one_frame();

                    tools.dbg('anp_gui[draw_run_continuous_animation]:\tdrawing %d :-)\n',this.a_time_ii);
                    pause(1/10);
                end
                tools.dbg('anp_gui[draw_run_continuous_animation]:\tStopping animation\n');

                this.s_draw_busy = false;
            catch err
                if ~strcmp(err.message,'Invalid or deleted object.')
                    rethrow(err);
                end
            end
        end
        
        function [] = draw_one_frame(this)
            this.draw_update_trails();
            this.draw_update_trail_head_arrows();
            this.draw_update_plot_titles();
            this.draw_update_textboxes();
        end
        
        function [] = draw_update_trails(this)
            % update the trail plot data
            set(this.h_z_plot_trail,'XData',real(this.d_z_values_truncated(this.d_t_trails{this.a_time_ii})),'YData',imag(this.d_z_values_truncated(this.d_t_trails{this.a_time_ii})));
            set(this.h_w_plot_trail,'XData',real(this.d_w_values_truncated(this.d_t_trails{this.a_time_ii})),'YData',imag(this.d_w_values_truncated(this.d_t_trails{this.a_time_ii})));        
        end
        
        function [] = draw_update_trail_head_arrows(this)
            % draw the (arrow-)head of the trail
            delete(this.h_z_plot_arrow); delete(this.h_w_plot_arrow);
            
            current_values_index =  this.a_time_ii * this.p_oversampling_factor;
            prev_values_index =     tools.iterator_modulo(this.a_time_ii * this.p_oversampling_factor - 1,this.p_n_time_steps * this.p_oversampling_factor);
            z_values_head_prev =    this.d_z_values(prev_values_index);
            w_values_head_prev =    this.d_w_values(prev_values_index);
            z_values_head =         this.d_z_values(current_values_index);
            w_values_head =         this.d_w_values(current_values_index);
            
            in_phi =                angle(z_values_head - z_values_head_prev);
            out_phi =               angle(w_values_head - w_values_head_prev);
            
            this.h_z_plot_arrow =       this.drawArrow(this.h_sub1,[real(this.d_z_values_truncated(this.a_time_ii * this.p_oversampling_factor)),imag(this.d_z_values_truncated(this.a_time_ii * this.p_oversampling_factor))],in_phi,this.p_z_arrow_length);
            this.h_w_plot_arrow =       this.drawArrow(this.h_sub2,[real(this.d_w_values_truncated(this.a_time_ii * this.p_oversampling_factor)),imag(this.d_w_values_truncated(this.a_time_ii * this.p_oversampling_factor))],out_phi,this.p_w_arrow_length);
            
        end
        
        function [] = draw_update_plot_titles(this)
            current_z_value = this.d_z_values(this.a_time_ii * this.p_oversampling_factor);
            current_w_value = this.d_w_values(this.a_time_ii * this.p_oversampling_factor);
            
            title(this.h_sub1,['Current value of D-Curve: ',num2str(current_z_value,'%.1f'),': M = ',num2str(abs(current_z_value),'%.2f'),' p = ',num2str(rad2deg(angle(current_z_value)),'%.2f'),'°']);
            title(this.h_sub2,['Nyquist: G(',num2str(current_z_value,'%.1f'),') = ',num2str(current_w_value,'%.1f'),': M = ',num2str(abs(current_w_value),'%.2f'),' p = ',num2str(rad2deg(angle(current_w_value)),'%.2f'),'°']);
        end
        
        % update the zero- and pole contribution and the cumulative values
        function [] = draw_update_textboxes(this)
            res_magnitude =                     '(';
            res_phase =                         '';
            
            for z = 1:this.d_n_zeros
                z_contribution =                this.d_z_values(this.a_time_ii * this.p_oversampling_factor) - this.d_zeros(z);
                this.h_text_z_annot(z).String = ['Z',num2str(z),': (',num2str(this.d_z_values(this.a_time_ii * this.p_oversampling_factor),           '%.1f'),') - (',num2str(this.d_zeros(z),'%.1f'),'): M=',num2str(abs(z_contribution),'%.1f'),' p=',num2str(rad2deg(angle(z_contribution)),'%.1f'),'°'];
                res_magnitude =                 [res_magnitude, '*',  num2str(abs(z_contribution),           '%.2f')];
                res_phase =                     [res_phase,     '+',  num2str(rad2deg(angle(z_contribution)),'%.2f')];
            end
            
            res_magnitude =                     [res_magnitude, ')/('];
            
            for p = 1:this.d_n_poles
                p_contribution =                this.d_z_values(this.a_time_ii * this.p_oversampling_factor) - this.d_poles(p);
                this.h_text_p_annot(p).String = ['P',num2str(p),': (',num2str(this.d_z_values(this.a_time_ii * this.p_oversampling_factor),           '%.1f'),') - (',num2str(this.d_poles(p),'%.1f'),'): M=',num2str(abs(p_contribution),'%.1f'),' p=',num2str(rad2deg(angle(p_contribution)),'%.1f'),'°'];            
                res_magnitude =                 [res_magnitude, '*',  num2str(abs(p_contribution),           '%.2f')];
                res_phase =                     [res_phase,     '-',  num2str(rad2deg(angle(p_contribution)),'%.2f')];
            end
            
            res_magnitude =                     [res_magnitude,') = ',num2str(abs(this.d_w_values(this.a_time_ii * this.p_oversampling_factor)),      '%.3f')];
            res_phase =                         [res_phase,    ') = ',num2str(rad2deg(angle(this.d_w_values(this.a_time_ii * this.p_oversampling_factor))),'%.3f')];

            this.h_text_res_annot(1).String =   ['Magnitude: ',  res_magnitude];
            this.h_text_res_annot(2).String =   ['Phase:       ',res_phase];
        end
        
        function [] = draw_update_limits_and_plots(this)
            subplot(this.h_sub1);
            this.p_z_xlim = xlim;
            this.p_z_ylim = ylim;

            this.calc_truncated_z_values();
            this.draw_update_full_z_plot();

            subplot(this.h_sub2);
            this.p_w_xlim = xlim;
            this.p_w_ylim = ylim;
            
            tools.dbg('anp_gui[draw_update_limits_and_plots]:\tz_xlim=[%.2f,%.2f], z_ylim =[%.2f,%.2f], w_xlim=[%.2f,%.2f], w_ylim=[%.2f,%.2f]\n',this.p_z_xlim(1),this.p_z_xlim(2),this.p_z_ylim(1),this.p_z_ylim(2),this.p_w_xlim(1),this.p_w_xlim(2),this.p_w_ylim(1),this.p_w_ylim(2));
            
            this.calc_truncated_w_values();
            this.draw_update_full_w_plot();
        end
        
        function [] = draw_update_full_z_plot(this)
            set(this.h_z_plot_full,'XData',real(this.d_z_values_truncated),'YData',imag(this.d_z_values_truncated));
        end
        
        function [] = draw_update_full_w_plot(this)
            set(this.h_w_plot_full,'XData',real(this.d_w_values_truncated),'YData',imag(this.d_w_values_truncated));
        end
        
        
        % -----------------------------------------------------------------
        % Callback methods
        % -----------------------------------------------------------------
        
        function [] = cb_step(this,~,~,dir)
            this.s_draw_busy = true;
            tools.dbg('anp_gui[cb_step]:\tStepping button.\n');
            
            this.a_time_ii = tools.iterator_modulo(this.a_time_ii + dir,this.p_n_time_steps);
            this.draw_one_frame();
            this.s_draw_busy = false;
        end
        
        
        function [] = cb_run(this,src,~,dir) % ignored parameters are src,evt
            if this.s_draw_busy && (dir == this.a_direction)
                return;
            elseif this.s_draw_busy && (dir ~= this.a_direction)
                this.a_direction = dir;
                for ii = 1:length(this.ui_run_switches)
                    if ~isequal(src,this.ui_run_switches(ii))
                        this.ui_run_switches(ii).State = 'off';
                    end
                end
            else
                tools.dbg('anp_gui[cb_run]:\tStart requested.\n');
                if this.s_data_ready
                    this.ui_sw_pause.State =    'off';
                    this.ui_btn_prev.Enable =   'off';
                    this.ui_btn_next.Enable =   'off';
                    for ii = 1:length(this.ui_run_switches)
                        if ~isequal(src,this.ui_run_switches(ii))
                            this.ui_run_switches(ii).State = 'off';
                        end
                    end
                    this.a_direction =          dir;
                    this.s_draw_allowed =       true;
                    this.draw_run_continuous_animation();
                else
                    warning('Data not ready!');
                    src.State =  'off';
                end
            end
        end
        
        function [] = cb_run_switch_off(this,src,~,dir)
            if this.s_draw_allowed == true && (dir == this.a_direction)
                src.State = 'on';
            end
        end
        
        function [] = cb_pause(this,src,~)
            tools.dbg('anp_gui[cb_pause]:\tPause requested.\n');
            this.s_draw_allowed =   false;
            src.State =             'on';
            for ii = 1:length(this.ui_run_switches)
                this.ui_run_switches(ii).State = 'off';
            end
            this.ui_btn_prev.Enable =   'on';
            this.ui_btn_next.Enable =   'on';
        end
        
        function cb_after_zoom_or_pan(this,~,~)
            if ~this.s_draw_busy
                
                this.draw_update_limits_and_plots();
                this.calc_plot_z_arrow_length();
                this.calc_plot_w_arrow_length();
                this.draw_one_frame();
            else
                this.s_check_limits =   true;
            end
        end
        
        % -----------------------------------------------------------------
        % utility functions
        % -----------------------------------------------------------------
        function arrowHandle = drawArrow(~,parent,x0,phi,l)
            r = l*[cos(phi),sin(phi)];
            
            arrowHandle = annotation('Arrow');
            set(arrowHandle,'parent',parent,'position',[x0-r,r]);
        end
        
        function h_arrow = draw_text_arrow(~,x0,phi,l,text,color)
            r = l*[cos(phi),sin(phi)];

            h_arrow = annotation('TextArrow');
            set(h_arrow,'parent',gca,'position',[x0-r,r],'String',text,'Color',color);
        end

        % this function clips the input values to the maximum value allowed inside
        % the plot
        function values_truncated = trunc(~,values,xlim,ylim)
            xlim = anp_stretch_centered(xlim,0.97);
            ylim = anp_stretch_centered(ylim,0.97);
            values_truncated = max(xlim(1), min(xlim(2), real(values))) ...
                          + 1i*max(ylim(1), min(ylim(2), imag(values)));
            %values_truncated = max(xlim(1), min(xlim(2), real(values))) ...
            %              + 1i*max(ylim(1), min(ylim(2), imag(values)));
        end
    end
end