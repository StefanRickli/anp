% -------------------------------------------------------------------------
%   Animated Nyquist Plot
%   *********************
%
%   Author: Stefan Rickli
%
%   Version 2.0, 2016-11-05
%   
%
%   What does it do?:
%   -----------------
%
%   This script calculates for a given set of zeros and poles a suitable
%   D-curve, plots it together with the zeros and poles and draws the
%   corresponding nyquist curve in another plot.
%   For every time step one is shown the contribution of each zero and pole
%   to the magnitude and phase of the nyquist curve's values.
%   After calculation, the animation is loaded into an interactive player
%   where one can further investigate each step.
%
%   known issues:
%   -------------
%
%   - axis scale sometimes disappears with certain combinations of
%     zeros/poles
%
% -------------------------------------------------------------------------

function [] = animated_nyquist_plot()
    clear; close all;
    
    %% debug
    t_min = 0;
    t_max = 1;
        
    %% user input: transfer function
    g.zeros = [-3];        % zeros of transfer function
    g.poles = [-1+1i,-1-1i];     % poles of transfer function
    
    %% initialization
    phi = 5;                                    % [°], roundoff-parameter of the D-curve, usually not necessary to change
    R = 2*max([abs(g.zeros) abs(g.poles)]);     % radius of the D-curve. auto value is such that one can oftentimes see what happens in origin of nyquist curve
    g.tf = zerosPoles2tfHandle(g.zeros, g.poles);
    animationParams = struct('duration', 0, 'fps', 0, 'trail_length', 0, 'x0_in', [], 'x0_out', [], 'dims_in', [], 'dims_out', [], 'N_frames', 0, 'N_trail', 0);
    
    
    %% user input: animation and window parameters
    animationParams.duration = 20;              % length of animation in movie player, affects spatial resolution, [s]
    animationParams.fps = 18;                   % speed of animation in movie player, affects spatial resolution, [frames/s]
    animationParams.trail_length = 0.15;        % how long should the trail be?, [0<trail_length<1]
    
    animationParams.in_auto_size = true;        % values below are only active if auto = false
    animationParams.in_x0 = [0;0];              % center of left plot
    animationParams.in_dims = [2.01*R;2.01*R];  % width and heigth of left plot
    
    animationParams.out_auto_size = true;       % values below are only active if auto = false
    animationParams.out_x0 = [0.2;0];           % center of right plot
    animationParams.out_dims = [2;2];           % width and heigth of right plot
    
    animationParams.plotSize = 500;             % pixel
    animationParams.border = 20;                % pixel   
    
    %% preparations
    % prepare frame buffer
    animationParams.N_frames = animationParams.duration * animationParams.fps;
    animationParams.N_trail = fix(animationParams.N_frames * animationParams.trail_length);
    
    % set up parametrisation according to duration and fps
    t_step = (t_max-t_min)/(animationParams.N_frames-1);
    t = t_min:t_step:t_max;
    
    % prepare t-intervals that are plotted on each frame, considering the
    % set number of trail-values
    t_indexes = cell(animationParams.N_frames,1);
    for ii = 1:animationParams.N_frames
        t_indexes{ii} = max(1,ii-animationParams.N_trail):ii;
    end
    
    %% calculate the function values
    in = getInputFunctionHandle();
    out = @(t,R,phi) g.tf(in(t,R,phi));
    
    in_values = in(t,R,phi);
    out_values = out(t,R,phi);
    
    %% calculation of the plot(s)
    frames = drawFunctions(in_values, out_values, g.zeros, g.poles, R, t_indexes, animationParams);
    
    %% display the movie player
    implay(frames,animationParams.fps);
    set(findall(0,'tag','spcui_scope_framework'),'position',[150 150 length(frames(1).cdata(1,:,1))+20 length(frames(1).cdata(:,1,1))+30]);
end

%% drawing function, returns a movie frame struct
function [frames] = drawFunctions(in_values, out_values, zeros, poles, R, t_indexes, animationParams)
    % remember how many zeros and poles we have
    zeros_N = length(zeros);
    poles_N = length(poles);

    % try to find optimal axis limits for the two plots, or let the user
    % decide if he chose to set either 'auto_size' to false
    [in_axis_xlim,in_axis_ylim,out_axis_xlim,out_axis_ylim] = autoZoom([zeros,poles],R,out_values);
    
    if ~animationParams.in_auto_size
       in_axis_xlim = [(animationParams.in_x0(1) - animationParams.in_dims(1)/2),(animationParams.in_x0(1) + animationParams.in_dims(1)/2)];
       in_axis_ylim = [(animationParams.in_x0(2) - animationParams.in_dims(2)/2),(animationParams.in_x0(2) + animationParams.in_dims(2)/2)];
    end
    
    if ~animationParams.out_auto_size
       out_axis_xlim = [(animationParams.out_x0(1) - animationParams.out_dims(1)/2),(animationParams.out_x0(1) + animationParams.out_dims(1)/2)];
       out_axis_ylim = [(animationParams.out_x0(2) - animationParams.out_dims(2)/2),(animationParams.out_x0(2) + animationParams.out_dims(2)/2)];
    end
    
    % derive the arrow length relative to subplot axis dimensions
    in_axis_width_1 = in_axis_xlim(2) - in_axis_xlim(1);        % [1]
    in_axis_height_1 = in_axis_ylim(2) - in_axis_ylim(1);       % [1]
    out_axis_width_1 = out_axis_xlim(2) - out_axis_xlim(1);     % [1]
    out_axis_height_1 = out_axis_ylim(2) - out_axis_ylim(1);    % [1]

    in_arrow_length = 0.04 * sqrt(in_axis_width_1^2 + in_axis_height_1^2);      % [1]
    out_arrow_length = 0.04 * sqrt(out_axis_width_1^2 + out_axis_height_1^2);   % [1]
    
    % preallocate space for frames to be stored
    frames(animationParams.N_frames) = struct('cdata',[],'colormap',[]);    
    
    % create the figure with its two subplots and remember their handles
    fig = figure;
    sub1 = subplot(1,2,1);
    hold on;
    sub2 = subplot(1,2,2);
    hold on;
    
    % set the visual params of the subplots
    subplot(sub1);
    axis equal; % for 1:1 aspect ratio
    xlim(in_axis_xlim), ylim(in_axis_ylim);
    ax = gca; ax.XAxisLocation = 'origin'; ax.YAxisLocation = 'origin'; % move the axes into the origin
    subplot(sub2);
    axis equal; % for 1:1 aspect ratio
    xlim(out_axis_xlim), ylim(out_axis_ylim);
    ax = gca; ax.XAxisLocation = 'origin'; ax.YAxisLocation = 'origin'; % move the axes into the origin
    
    % calculate the correct positions (relative to interior of the figure)
    fig_plot_height = 2*animationParams.border + animationParams.plotSize;              % plot + border above and below, pixel
    fig_annotation_textbox_height = 14;                                                 % one box, pixel
    fig_annotation_height_sum = (max(zeros_N,poles_N)+2)*fig_annotation_textbox_height; % cumulative, pixel
    fig_height = fig_plot_height + fig_annotation_height_sum;                           % pixel
    fig_width = 3*animationParams.border + 2*animationParams.plotSize;                                  % pixel
    
    % fig.Position expects pixels as unit
    fig.Position = [100 100 fig_width fig_height];
    
    % as subplot.Postion expects fractions of the inside of the figure, we
    % recalculate the pixel values
    fig_border_horizontal_frac = animationParams.border/fig_width;                      % fraction
    fig_border_vertical_frac = animationParams.border/fig_height;                       % fraction
    fig_plot_width_frac = (fig_width - 3*animationParams.border)/(2*fig_width);         % fraction
    fig_plot_height_frac = 1-(2*animationParams.border + fig_annotation_height_sum)/fig_height; % fraction
    
    % subX.Position expects fractions of the inside of the figure as unit
    sub1.Position = [fig_border_horizontal_frac, (1 - fig_border_vertical_frac - fig_plot_height_frac), fig_plot_width_frac, fig_plot_height_frac];
    sub2.Position = [(2*fig_border_horizontal_frac + fig_plot_width_frac),(1 - fig_border_vertical_frac - fig_plot_height_frac), fig_plot_width_frac, fig_plot_height_frac];
    
    % draw the zeros and poles in the left (input function) subplot
    subplot(sub1);
    for z = 1:zeros_N
        zero = zeros(z);
        zero_trunc = trunc(zero, in_axis_xlim, in_axis_ylim);
        if zero ~= zero_trunc
            % The zero in question being not equal to its truncated value
            % means that it lies outside the current subplot limits.
            % Instead of drawing just a circle at the border, plot an
            % arrow, indicating that there's an outlier.
            drawTextArrow([real(zero_trunc),imag(zero_trunc)],angle(zero),in_arrow_length,' o',[95 158 160]/255);
        else
            scatter(real(zero),imag(zero),60,'o','MarkerEdgeColor',[70 130 180]/255,'LineWidth',1.5);
        end
    end

    for p = 1:poles_N
        pole = poles(p);
        pole_trunc = trunc(pole, in_axis_xlim, in_axis_ylim);
        if pole ~= pole_trunc
            % same conclusion as above with the zeros
            drawTextArrow([real(pole_trunc),imag(pole_trunc)],angle(pole),in_arrow_length,' x',[255 140 0]/255);
        else
            scatter(real(pole),imag(pole),60,'x','MarkerEdgeColor',[255 140 0]/255,'LineWidth',1.5);
        end
    end
    
    % calculate the truncated values of the function (according to window
    % params)
    in_values_truncated = trunc(in_values, in_axis_xlim, in_axis_ylim);
    out_values_truncated = trunc(out_values, out_axis_xlim, out_axis_ylim);

    % plot the full input- and output curves with some transparency
    in_plot_full = plot(sub1,in_values_truncated,'Color',[0.05 0.4970 0.7410]); in_plot_full.Color(4) = 0.3;
    out_plot_full = plot(sub2,out_values_truncated,'Color',[0.05 0.4970 0.7410]); out_plot_full.Color(4) = 0.3;
    
    % prepare the curves' trail plots and remember their handle for later
    % use
    in_plot_trail = plot(sub1,0,0); set(in_plot_trail,'linewidth',2,'Color',[255 215 0]/255);
    out_plot_trail = plot(sub2,0,0); set(out_plot_trail,'linewidth',2,'Color',[255 215 0]/255);

    % prepare the arrows annotations that mark the value of the head of the
    % trail at the current frame and remember their handle for later use
    subplot(sub1); in_plot_arrow = annotation('line',[0 0],[1 0]);
    subplot(sub2); out_plot_arrow = annotation('line',[0 0],[1 0]);
    
    % We need to know the value of the head of the trail at the previous
    % frame in order to determine the arrow's direction. In the case of the
    % start direction, we need the value at the second to last frame
    % because if we took the last frame then we'd use f(t=1) = f(t=0).
    in_values_head_prev = in_values(animationParams.N_frames - 1);
    out_values_head_prev = out_values(animationParams.N_frames - 1);
    
    
    % prepare the annotations below the plots and remember their handles
    
    % again we need to know relative (fractions) positions inside the figure
    fig_annotation_start_frac = (fig_annotation_height_sum + 10)/fig_height;    % fraction
    fig_annotation_textbox_frac = fig_annotation_textbox_height/fig_height;     % fraction
    
    figure(fig);
    % zero and pole contributions
    for ii = 1:zeros_N
        z_annotations(ii) = annotation('TextBox',[fig_border_horizontal_frac, (fig_annotation_start_frac - fig_annotation_textbox_frac*(ii-1)), fig_plot_width_frac, fig_annotation_textbox_frac],'String',['zero ',num2str(ii)],'LineStyle','none','FontSize',9);
    end
    for ii = 1:poles_N
        p_annotations(ii) = annotation('TextBox',[0.5, (fig_annotation_start_frac - fig_annotation_textbox_frac*(ii-1)), fig_plot_width_frac, fig_annotation_textbox_frac],'String',['pole ',num2str(ii)],'LineStyle','none','FontSize',9);
    end
    
    % cumulative calculations
    res_annotation1 = annotation('TextBox',[fig_border_horizontal_frac max(0,(fig_annotation_start_frac - fig_annotation_textbox_frac*(max(zeros_N,poles_N)+1))) fig_plot_width_frac fig_annotation_textbox_frac],'String','resultline 1','LineStyle','none','FontSize',9);
    res_annotation2 = annotation('TextBox',[fig_border_horizontal_frac max(0,(fig_annotation_start_frac - fig_annotation_textbox_frac*(max(zeros_N,poles_N)+2))) fig_plot_width_frac fig_annotation_textbox_frac],'String','resultline 2','LineStyle','none','FontSize',9);
    
    % main loop where we update the plots and annotations and store the
    % figure in a movie frame after each iteration
    for ii = 1:animationParams.N_frames        
        % update the trail plot data
        set(in_plot_trail,'XData',real(in_values_truncated(t_indexes{ii})),'YData',imag(in_values_truncated(t_indexes{ii})));
        set(out_plot_trail,'XData',real(out_values_truncated(t_indexes{ii})),'YData',imag(out_values_truncated(t_indexes{ii})));        
        
        % draw the (arrow-)head of the trail
        delete(in_plot_arrow); delete(out_plot_arrow);
        in_values_head = in_values(ii);
        out_values_head = out_values(ii);
        in_phi = angle(in_values_head - in_values_head_prev);
        out_phi = angle(out_values_head - out_values_head_prev);
        subplot(sub1); in_plot_arrow = drawArrow([real(in_values_truncated(ii)),imag(in_values_truncated(ii))],in_phi,in_arrow_length);
        subplot(sub2); out_plot_arrow = drawArrow([real(out_values_truncated(ii)),imag(out_values_truncated(ii))],out_phi,out_arrow_length);
        
        % update the title
        title(sub1,['value at arrow = ',num2str(in_values_head,'%.1f'),': M = ',num2str(abs(in_values_head),'%.2f'),' p = ',num2str(rad2deg(angle(in_values_head)),'%.2f'),'°']);
        title(sub2,['value at arrow = ',num2str(out_values_head,'%.1f'),': M = ',num2str(abs(out_values_head),'%.2f'),' p = ',num2str(rad2deg(angle(out_values_head)),'%.2f'),'°']);
        
        % update the zero- and pole contribution and the cumulative values
        res_magnitude = '(';
        res_phase = '';
        for z = 1:zeros_N
            z_contribution = in_values(ii) - zeros(z);
            z_annotations(z).String = ['Z',num2str(z),': (',num2str(in_values(ii),'%.1f'),') - (',num2str(zeros(z),'%.1f'),'): M=',num2str(abs(z_contribution),'%.1f'),' p=',num2str(rad2deg(angle(z_contribution)),'%.1f'),'°'];
            res_magnitude = [res_magnitude,'*',num2str(abs(z_contribution),'%.2f')];
            res_phase = [res_phase,'+',num2str(rad2deg(angle(z_contribution)),'%.2f')];
        end
        res_magnitude = [res_magnitude,')/('];
        for p = 1:poles_N
            p_contribution = in_values(ii) - poles(p);
            p_annotations(p).String = ['P',num2str(p),': (',num2str(in_values(ii),'%.1f'),') - (',num2str(poles(p),'%.1f'),'): M=',num2str(abs(p_contribution),'%.1f'),' p=',num2str(rad2deg(angle(p_contribution)),'%.1f'),'°'];            
            res_magnitude = [res_magnitude,'*',num2str(abs(p_contribution),'%.2f')];
            res_phase = [res_phase,'-',num2str(rad2deg(angle(p_contribution)),'%.2f')];
        end
        res_magnitude = [res_magnitude,') = ',num2str(abs(out_values(ii)),'%.3f')];
        res_phase = [res_phase,') = ',num2str(rad2deg(angle(out_values(ii))),'%.3f')];
        
        res_annotation1.String = ['Magn.: ',res_magnitude];
        res_annotation2.String = ['Phase: ',res_phase];
        
        % save the current frame
        frames(ii) = getframe(gcf);
        
        in_values_head_prev = in_values_head;
        out_values_head_prev = out_values_head;        
    end
    close(fig);
end

% The purpose of this function is to provide a means of detecting clumping
% of zeros and poles. If we have a bunch of zeros and poles that are fairly
% close to each other with one outlier, forget about the outlier and zoom
% in on the group.
function [in_axis_xlim,in_axis_ylim,out_axis_xlim,out_axis_ylim] = autoZoom(zp, R, out_values)
    % find the 'center of gravity' of the zeros and poles
    locus_m = mean(zp);
    
    % where are 66% of the points?
    locus_std_real = std(real(zp));
    locus_std_imag = std(imag(zp));
    
    % get indexes of those points that lie within the ellipse given by the stddev in x- and y-direction
    locus_inner =  find(abs(zp-locus_m) <= sqrt(locus_std_real^2+locus_std_imag^2));
    
    % how much area do these points cover in the zoomed-out plot (when the D-curve fits into it)?
    locus_inner_dims = [max(real(zp(locus_inner)))-min(real(zp(locus_inner)));max(imag(zp(locus_inner)))-min(imag(zp(locus_inner)))];
    locus_inner_area = locus_inner_dims(1)*locus_inner_dims(2);
    locus_fract = locus_inner_area/(2.01*R)^2;
    
    % if that area is below 0.2% then we zoom into the bulk of the points,
    % leaving out outliers
    if locus_fract > 0.002
        % all good, fit the whole D-curve and the zp-points into the plot
        in_axis_xlim = 1.05*R*[-1,1];
        in_axis_ylim = 1.05*R*[-1,1];
    else
        % include the origin in the box around the clumped points and then
        % make the plot size 8 times (rule of thumb) larger
        in_axis_xlim = resize_limits([min([real(zp(locus_inner)),0]),max([real(zp(locus_inner)),0])],8);
        in_axis_ylim = resize_limits([min([imag(zp(locus_inner)),0]),max([imag(zp(locus_inner)),0])],8);
        
        width = in_axis_xlim(2) - in_axis_xlim(1);
        height = in_axis_ylim(2) - in_axis_ylim(1);
        
        % make the plot rectangular (make the smaller axis as large as the
        % bigger one)
        if width >= height
            y0 = mean(in_axis_ylim);
            in_axis_ylim = [(y0-width/2),(y0+width/2)];
            height = width;
        else
            x0 = mean(in_axis_xlim);
            in_axis_xlim = [(x0-height/2),(x0+height/2)];
            width = height;
        end
    end
    
    % for the right subplot it's fairly easy. Have the Plot such that the
    % whole nyquist curve fits into it.
    out_axis_xlim = resize_limits([min(real(out_values)),max(real(out_values))],1.1);
    out_axis_ylim = resize_limits([min(imag(out_values)),max(imag(out_values))],1.1);
    
    width = out_axis_xlim(2) - out_axis_xlim(1);
    height = out_axis_ylim(2) - out_axis_ylim(1);

    % make this plot also rectangular
    if width >= height
        y0 = mean(out_axis_ylim);
        out_axis_ylim = [(y0-width/2),(y0+width/2)];
        height = width;
    else
        x0 = mean(out_axis_xlim);
        out_axis_xlim = [(x0-height/2),(x0+height/2)];
        width = height;
    end

end

function arrowHandle = drawArrow(x0, phi, l)
    r = l*[cos(phi),sin(phi)];
    
    arrowHandle = annotation('Arrow');
    set(arrowHandle,'parent',gca,'position',[x0-r,r]);
end

function arrowHandle = drawTextArrow(x0, phi, l, text, color)
    l = 0.2;                    % length of arrow
    r = l*[cos(phi),sin(phi)];
    
    arrowHandle = annotation('TextArrow');
    set(arrowHandle,'parent',gca,'position',[x0-r,r],'String',text,'Color',color);
end

%% input function shapes
% t is in the following range 0<=t<=1.
% the shape should be continuous and a closed curve (i.e. shape(0)=shape(1))
function inputFunctionHandle = getInputFunctionHandle()
    %inputFunctionHandle = @(t) unitCircle(t);
    %inputFunctionHandle = @(t) dCurve(t);
    inputFunctionHandle = @(t,R,phi) roundedCurve(t,R,phi);
end

% unused, test function
function z = unitCircle(x)
    z = exp(2*pi*i.*x);
end

% this D-shaped curve has small arcs at the edges to smooth them and
% make the derivative of the curve continuous
% The shape is basically split into 5 parts to which the indexes correspond
% The two parameters one can play with are phi, which determines the degree
% of roundoff (0<phi<180) and the radius of the big arc R.
function z = roundedCurve(t,R,phi)
    phi = deg2rad(phi);
    
    % calculate values of the small round-off circles
    r = R*sin(phi)/(1+sin(phi));
    m2 = r + i*sqrt(R^2-2*R*r);
    m4 = r - i*sqrt(R^2-2*R*r);
    
    % need to know how long the parts are
    d1 = sqrt(R^2-2*R*r);
    d2 = r*(pi/2 + phi);
    d3 = R*(pi - 2*phi);
    d4 = d2;
    d5 = d1;
    
    % length of the round part (small arcs + big arc)
    d_arcs = d2 + d3 + d4;
    
    % calculate the t-values at which two adjacent regions touch
    % the only arbitrary choice here is t12. together with t45 it
    % determines the fraction of time we need for the y-axis part of the
    % D-curve.
    % the v-values are velocity-factors
    t12 = 2/5; 
    t45 = 1-t12;
    v_ave = d_arcs / (t45 - t12);
    
    t23 = t12 + d2/v_ave;
    t34 = 1-t23;
    
    v11 = 3; % again the only arbitrary choice amongst the velocity-factors. determines the speed of the animation near the origin.
    v12 = (d1-v11*t12)/t12^4;
    v51 = v11;
    v52 = v12;
    
    z = zeros(length(t),1);
    for ii = 1:length(t)
        if t(ii) <= t12
            z(ii) = (v11*t(ii) + v12*t(ii)^4)*i; % part 1
        elseif t(ii) > t12 && t(ii) < t23
            z(ii) = m2 + r*exp(i*(pi - (t(ii)-t12)*(pi/2 + phi)/(t23-t12))); % part 2
        elseif t(ii) > t23 && t(ii) < t34
            z(ii) = R*exp(i*(pi/2 - phi - (t(ii)-t23)*(pi-2*phi)/(t34-t23))); % part 3
        elseif t(ii) > t34 && t(ii) < t45
            z(ii) = m4 + r*exp(i*(phi - pi/2 - (t(ii)-t34)*(pi/2 + phi)/(t45-t34))); % part 4
        else
            z(ii) = -(v51*(1-t(ii)) + v52*(1-t(ii))^4)*i; % part 5
        end
    end
end

% unused, this D-curve has sharp edges which cause the arrow in the plot to jump!
function z = dCurve(x)
    r = 10;
    part = 0.3;
    v = 90;
    t = 3;
    c = (r-t*part)/exp(v*part);
    
    z = zeros(length(x),1);
    for ii = 1:length(x)
        if x(ii) <= part
            z(ii) = (t*x(ii)+c*exp(v*x(ii)))*i;
        elseif x(ii) > part && x(ii) < 1-part
            z(ii) = r*exp(i* (pi/2 * (1-(x(ii)-part)/(0.5-part))));
        else
            z(ii) = -(t*(1-x(ii))+c*exp(v*(1-x(ii))))*i;
        end
    end
end

%% transfer function preparation
function g = zerosPoles2tfHandle(zeros,poles)
    g = @(t) polyval(poly(zeros),t)./polyval(poly(poles),t);
end

%% various stuff
% this function clips the input values to the maximum value allowed inside
% the plot
function z = trunc(x,xlim,ylim)
    z = max(xlim(1), min(xlim(2), real(x))) ...
        + i*max(ylim(1), min(ylim(2), imag(x)));
end

% stretches a given interval with a given factor about its center
function result = resize_limits(in, factor)
    size_in = in(2) - in(1);
    result = [(in(1)-(factor-1)*size_in/2),(in(2)+(factor-1)*size_in/2)];
end