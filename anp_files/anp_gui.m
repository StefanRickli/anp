% as in p.41
% last line: 516
classdef anp_gui < handle
    properties(SetAccess = private)
        % g: general properties
        g_uid;
        g_legacy            string;
        
        % s: state
        s_draw_busy         logical;
        
        % h: handles
        h_anp_tf_processor
        
        h_fig               matlab.ui.Figure;
        h_sub1              matlab.graphics.axis.Axes;
        h_sub2              matlab.graphics.axis.Axes;
        
        h_text_z_annot      matlab.graphics.shape.TextBox
        h_text_p_annot      matlab.graphics.shape.TextBox
        h_text_res_annot    matlab.graphics.shape.TextBox
        
        h_z_plot_full       matlab.graphics.chart.primitive.Line
        h_z_plot_trail      matlab.graphics.chart.primitive.Line
        h_z_plot_arrow      matlab.graphics.shape.Arrow
        h_z_pz_objcts     % array of handles with possibly mixed type
        
        h_w_plot_full       matlab.graphics.chart.primitive.Line
        h_w_plot_trail      matlab.graphics.chart.primitive.Line
        h_w_plot_arrow      matlab.graphics.shape.Arrow
        
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
        p_time
        
        p_z_x0              double
        p_z_dims            double
        p_z_auto_lims       logical
        p_z_xlim            double
        p_z_ylim            double
        p_z_width           double
        p_z_height          double
        p_z_arrow_length    double
        
        p_w_x0              double
        p_w_dims            double
        p_w_auto_lims       logical
        p_w_xlim            double
        p_w_ylim            double
        p_w_width           double
        p_w_height          double
        p_w_arrow_length    double
        
        % d: data
        d_valid                 logical
        d_n_poles               double
        d_n_zeros               double
        d_poles                 double
        d_zeros                 double
        
        d_t_values              double
        d_t_oversampled         double
        
        d_z_values              double
        d_z_values_truncated    double
        d_z_values_head_prev    double
        
        d_w_values              double
        d_w_values_truncated    double
        d_w_values_head_prev    double
        
        d_radii
        
    end
    methods(Access = public)
        function this = anp_gui()
            this.g_uid =        java.util.UUID.randomUUID.toString;
            this.g_legacy =     anp_check_Matlab_version();
            this.h_fig =        figure('CloseRequestFcn',@this.on_figure_delete);
            this.h_sub1 =       subplot(1,2,1);
            hold on;
            this.h_sub2 =       subplot(1,2,2);
            hold on;
            fprintf('anp_gui\t\t\t\t%s: Instance created.\n',this.g_uid);
        end
        
        function [] = set_window_props(this,new_props)
            this.w_border =     new_props.border;
            this.w_plot_size =  new_props.plot_size;
            % TODO everything window and generally GUI related, NOT
            % plotting
        end
        
        function [] = set_z_plot_limits(this,new_props)
            % TODO check for proper array dimensions
            
            this.p_z_x0 =           new_props.z_x0;
            this.p_z_dims =         new_props.z_dims;
            this.p_z_auto_lims =    new_props.z_auto_lims;
        end
        
        function [] = set_w_plot_limits(this,new_props)
            % TODO check for proper array dimensions
            
            this.p_w_x0 =           new_props.w_x0;
            this.p_w_dims =         new_props.w_dims;
            this.p_w_auto_lims =    new_props.w_auto_lims;
        end
        
        function [] = set_radii(this,new_props)
            % TODO check for proper array dimensions
            
            this.d_radii =          new_props.radii;
        end
        
        function [] = set_poles_zeros(this,new_props)
            % TODO check for proper array dimensions
            
            this.d_poles =          new_props.poles;
            this.d_zeros =          new_props.zeros;
            
            this.d_n_poles =        length(new_props.poles);
            this.d_n_zeros =        length(new_props.zeros);
        end

        function [] = set_anp_tf_processor(this,h)
            this.h_anp_tf_processor = h;
        end
        
        function [] = init_visuals(this)
            this.fetch_data();
            this.calc_gui_positions();
            this.calc_plot_axis_limits();
            this.calc_plot_w_arrow_length();
            this.calc_plot_z_arrow_length();
            this.calc_truncated_z_values();
            this.calc_truncated_w_values();
            this.draw_init_gui_statics();
            this.draw_init_gui_text_objects();
            this.draw_init_plot_axes();
            this.draw_init_plot_arrows();
            this.draw_init_z_plot_poles_zeros();
        end
        
        function delete(this)
            fprintf('anp_gui\t\t\t\t%s: Deletion requested.\n',this.g_uid);
            if isvalid(this.h_fig)
                close(this.h_fig);
            end
        end
    end
    methods(Access = private)
        function on_figure_delete(this,src,~) % last argument is 'evt'
            % TODO remove all event listeners here!
            disp('Deleting figure window.');
            delete(src);
            this.delete();
        end
        
        function [] = fetch_data(this)
            assert(isvalid(this.h_anp_tf_processor));
            
            time_data = this.h_anp_tf_processor.get_time_points();
            this.d_t_values =       time_data.time_points;
            this.d_t_oversampled =  time_data.data_points;
            this.p_time =           time_data.time_props;
            
            function_data = this.h_anp_tf_processor.get_data();
            this.d_z_values =       function_data.z_values;
            this.d_w_values =       function_data.w_values;
        end
        
        function [] = draw_init_gui_statics(this)
            % TODO 
            this.h_fig.Position =   this.w_fig_position;
            
            this.h_sub1.Position =  this.w_sub1_position;
            this.h_sub2.Position =  this.w_sub2_position;
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
            this.h_z_plot_trail =           plot(sub1,0,0);
            set(this.h_z_plot_trail,'Color',[255 215 0]/255,'linewidth',2);
            this.h_w_plot_trail =           plot(sub2,0,0);
            set(this.h_w_plot_trail,'Color',[255 215 0]/255,'linewidth',2);
        end
        
        function [] = draw_update_full_z_plot(this)
            set(this.h_z_plot_full,'XData',real(this.d_z_values_truncated),'YData',imag(this.d_z_values_truncated));
        end
        function [] = draw_update_full_w_plot(this)
            set(this.h_w_plot_full,'XData',real(this.d_w_values_truncated),'YData',imag(this.d_w_values_truncated));
        end
        
        function [] = draw_init_plot_arrows(this)
            % prepare the arrows annotations that mark the value of the head of the
            % trail at the current frame and remember their handle for later use
            subplot(this.h_sub1);
            this.h_z_plot_arrow =       annotation('Arrow',[0 0],[1 0]);

            subplot(this.h_sub2);
            this.h_w_plot_arrow =       annotation('Arrow',[0 0],[1 0]);
            
            % We need to know the value of the head of the trail at the previous
            % frame in order to determine the arrow's direction.
            this.d_z_values_head_prev = this.d_z_values(this.p_time.oversampling_factor-1);
            this.d_w_values_head_prev = this.d_w_values(this.p_time.oversampling_factor-1);
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
        
        function [] = draw_init_plot_axes(this)
            % TODO maybe separate into two separate functions for z and w
            % plot
            
            subplot(this.h_sub1);
            axis equal; % for 1:1 aspect ratio
            xlim(this.p_z_xlim), ylim(this.p_z_ylim);
            anp_plot_axes_origin(this.g_legacy);
            grid on;
            
            subplot(this.h_sub2);
            axis equal;
            xlim(this.p_w_xlim), ylim(this.p_w_ylim);
            anp_plot_axes_origin(this.g_legacy);
            grid on;
        end
        
        % draw the poles and zeros in the left (input function) subplot
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
        
        function h_arrow = draw_text_arrow(x0,phi,l,text,color)
            r = l*[cos(phi),sin(phi)];

            h_arrow = annotation('TextArrow');
            set(h_arrow,'parent',gca,'position',[x0-r,r],'String',text,'Color',color);
        end
        
        function [] = calc_plot_axis_limits(this)
            % TODO make sure that all props and data are set before calling
            % this method
            
            assert(~any([isempty(this.d_zeros),...
                         isempty(this.d_poles),...
                         isempty(this.d_radii),...
                         isempty(this.d_w_values),...
                         isempty(this.p_z_auto_lims),...
                         isempty(this.p_w_auto_lims)]));
            
            if this.p_z_auto_lims
                % try to find optimal axis limits for the two plots, or let the user
                % decide if he chose to set either 'auto_size' to false
                
                [this.p_z_xlim,this.p_z_ylim] = anp_plot_auto_zoom_z([this.d_zeros,this.d_poles],this.d_radii.R);
            else
               this.p_z_xlim = [(this.p_z_x0(1) - this.p_z_dims(1)/2),(this.p_z_x0(1) + this.p_z_dims(1)/2)];
               this.p_z_ylim = [(this.p_z_x0(2) - this.p_z_dims(2)/2),(this.p_z_x0(2) + this.p_z_dims(2)/2)];
            end
            
            if this.p_w_auto_lims
                [this.p_w_xlim,this.p_w_ylim] = anp_plot_auto_zoom_w(this.d_w_values);
            else
               this.p_w_xlim = [(this.p_w_x0(1) - this.p_w_dims(1)/2),(this.p_w_x0(1) + this.p_w_dims(1)/2)];
               this.p_w_ylim = [(this.p_w_x0(2) - this.p_w_dims(2)/2),(this.p_w_x0(2) + this.p_w_dims(2)/2)];
            end
        end
        
        function [] = calc_plot_z_arrow_length(this)
            in_axis_width = diff(this.p_z_xlim);        % [1]
            in_axis_height = diff(this.p_z_ylim);       % [1]
            this.p_z_arrow_length = 0.04 * sqrt(in_axis_width^2 + in_axis_height^2);      % [1]
        end
        function [] = calc_plot_w_arrow_length(this)
            % TODO handle very short arrow lengths better
            out_axis_width = diff(this.p_w_xlim);     % [1]
            out_axis_height = diff(this.p_w_ylim);    % [1]
            this.p_w_arrow_length = 0.04 * sqrt(out_axis_width^2 + out_axis_height^2);   % [1]
            if this.p_w_arrow_length < 0.05
                this.p_w_arrow_length = 0.0001;
            end
        end
        
        function [] = calc_truncated_z_values(this)
            this.d_z_values_truncated =  this.trunc(this.d_z_values,this.p_z_xlim,this.p_z_ylim);
        end
        function [] = calc_truncated_w_values(this)
            this.d_w_values_truncated =  this.trunc(this.d_w_values,this.p_w_xlim,this.p_w_ylim);
        end
        % this function clips the input values to the maximum value allowed inside
        % the plot
        function values_truncated = trunc(~,values,xlim,ylim)
            values_truncated = max(xlim(1), min(xlim(2), real(values))) ...
                          + 1i*max(ylim(1), min(ylim(2), imag(values)));
        end
    end
end