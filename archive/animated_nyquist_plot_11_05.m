% -------------------------------------------------------------------------
%   Animated Nyquist Plot
%   *********************
%
%   Author: Stefan Rickli
%
%   Version 1.0, 2016-11-05
%   
% -------------------------------------------------------------------------

function [] = animated_nyquist_plot()
    clear; close all;
    
    %% debug
    t_min = 0;
    t_max = 1;
    
    %% user input: transfer function
    g.zeros = [-0.5];        % zeros of transfer function
    g.poles = [-1,-1,-1];     % poles of transfer function
    
    %% initialization
    g.tf = getTransferFunction(g.zeros, g.poles);
    global animationParams;
    animationParams = struct('duration', 0, 'fps', 0, 'trail_length', 0, 'x0_in', [], 'x0_out', [], 'dims_in', [], 'dims_out', [], 'N_frames', 0, 'N_trail', 0);
    
    
    %% user input: animation and window parameters
    animationParams.duration = 20;      % [s]
    animationParams.fps = 20;           % [frames/s]
    animationParams.trail_length = 0.1; % [0<trail_length<1]
    
    animationParams.x0_in = [0;0];      % center of left plot
    animationParams.x0_out = [0.2;0];   % center of right plot
    animationParams.dims_in = [5;5];    % width and heigth of left plot
    animationParams.dims_out = [2;2];   % width and heigth of right plot   
    
    %% preparations
    % prepare frame buffer
    animationParams.N_frames = animationParams.duration * animationParams.fps;
    animationParams.N_trail = fix(animationParams.N_frames * animationParams.trail_length);
    
    % set up parametrisation according to duration and fps
    t_step = (t_max-t_min)/(animationParams.N_frames-1);
    t = t_min:t_step:t_max;
    
    % prepare t-intervals that are plotted on each frame, considering the
    % set number of trail-values
    t_vec = cell(animationParams.N_frames,1);
    for ii = 1:animationParams.N_frames
        t_vec{ii} = t(max(1,ii-animationParams.N_trail):ii);
    end
    
    t_indexes = cell(animationParams.N_frames,1);
    for ii = 1:animationParams.N_frames
        t_indexes{ii} = max(1,ii-animationParams.N_trail):ii;
    end
    
    %% calculate the function values
    in = getInputFunctionHandle();
    nyq = @(t) g.tf(in(t));
    
    in_values = in(t);
    nyq_values = nyq(t);
    
    %% drawings
    global windowLimits;
    
    % plot of the transfer function (nyquist plot)
    windowLimits.xmin = animationParams.x0_out(1) - animationParams.dims_out(1)/2;
    windowLimits.xmax = animationParams.x0_out(1) + animationParams.dims_out(1)/2;
    windowLimits.ymin = animationParams.x0_out(2) - animationParams.dims_out(2)/2;
    windowLimits.ymax = animationParams.x0_out(2) + animationParams.dims_out(2)/2;
    
    frames_nyquist = drawFunction(nyq_values, g.zeros, g.poles, t, t_indexes, animationParams.x0_out, animationParams.dims_out);
    
    % plot of the corresponding input function (D-curve)
    windowLimits.xmin = animationParams.x0_in(1) - animationParams.dims_in(1)/2;
    windowLimits.xmax = animationParams.x0_in(1) + animationParams.dims_in(1)/2;
    windowLimits.ymin = animationParams.x0_in(2) - animationParams.dims_in(2)/2;
    windowLimits.ymax = animationParams.x0_in(2) + animationParams.dims_in(2)/2;
    
    frames_input = drawFunction(in_values, g.zeros, g.poles, t, t_indexes, animationParams.x0_in, animationParams.dims_in);    
    
    %% combining the frames
    frames_combined(animationParams.N_frames) = struct('cdata',[],'colormap',[]);
    for ii = 1:animationParams.N_frames
        frames_combined(ii).cdata = zeros(size(frames_input(1).cdata)*[1,0,0;0,2,0;0,0,1],'uint8');
        frames_combined(ii).cdata(:,1:length(frames_input(1).cdata(1,:,1)),:) = frames_input(ii).cdata;
        frames_combined(ii).cdata(:,length(frames_input(1).cdata(1,:,1))+1:2*length(frames_input(1).cdata(1,:,1)),:) = frames_nyquist(ii).cdata;
    end
    
    %% display the combined frames
    figure;
    imshow(frames_combined(1).cdata); % only there to circumvent a display bug
    movie(frames_combined,100,animationParams.fps);
    
    
end

%% drawing functions
function [frames] = drawFunction(values, zeros, poles, t, t_indexes, x0, dims)
    
    % calculate the truncated values of the function (according to window
    % params)
    truncated_values = trunc(values);
    
    % preallocate space for frames to be stored
    global animationParams;
    frames(animationParams.N_frames) = struct('cdata',[],'colormap',[]);    
    
    fig = figure;%('Visible', 'off');
    hold on;
    x_params = [x0(1) - dims(1)/2, x0(1) + dims(1)/2];
    y_params = [x0(2) - dims(2)/2, x0(2) + dims(2)/2];
    
    % plot zeros and poles
    N_zeros = length(zeros);
    for z = 1:N_zeros;
        scatter(real(zeros(z)),imag(zeros(z)),'o');
    end
    
    N_poles = length(poles);
    for p = 1:N_poles;
        scatter(real(poles(p)),imag(poles(p)),'x');
    end
    
    % plot the full curve with some transparency
    p_full = plot(truncated_values);
    p_full.Color(4) = 0.3;
    
    % prepare intialized objects to work on during the for-loop
    p_trail = plot(0,0);
    set(p_trail,'linewidth',2);
    
    arrow = annotation('line',[0 0],[1 0]);
    
    head_value_prev = values(animationParams.N_frames - 1); % needed for arrow direction

    % set the visual params of the figure
    axis equal; % for 1:1 aspect ratio
    xlim(x_params), ylim(y_params);
    ax = gca; ax.XAxisLocation = 'origin'; ax.YAxisLocation = 'origin'; % move the axes into the origin
    
    for ii = 1:animationParams.N_frames        
        % plot the trail
        set(p_trail,'XData',real(truncated_values(t_indexes{ii})),'YData',imag(truncated_values(t_indexes{ii})));
        
        % draw the (arrow-)head of the animation
        delete(arrow);        
        head_value = values(ii);
        phi = angle(head_value - head_value_prev);
        arrow = drawArrow([real(truncated_values(ii)),imag(truncated_values(ii))],phi);
        
        title(['value at triangle = ' num2str(head_value)]);
        
        % save the current frame
        frames(ii) = getframe(gcf);
        
        head_value_prev = head_value;
    end
    close(fig);
end

function arrowHandle = drawArrow(x0, phi)
    l = 0.2;                    % length of arrow
    r = l*[cos(phi),sin(phi)];
    
    arrowHandle = annotation('Arrow');
    set(arrowHandle,'parent',gca,'position',[x0-r,r]);
end

%% shapes
% t is in the following range 0<=t<=1.
% the shape should be continuous and a closed curve (i.e. shape(0)=shape(1))
function inputFunctionHandle = getInputFunctionHandle()
    %inputFunctionHandle = @(t) unitCircle(t);
    %inputFunctionHandle = @(t) dCurve(t);
    inputFunctionHandle = @(t) roundedCurve(t);
end

% testfunction
function z = unitCircle(x)
    z = exp(2*pi*i.*x);
end

% this D-shaped curve has small arcs at the edges to smooth them and
% make the derivative of the curve continuous
% The shape is basically split into 5 parts to which the indexes correspond
% The two parameters one can play with are phi, which determines the degree
% of roundoff (0<phi<180) and the radius of the big arc R.
function z = roundedCurve(t)
    % round-off parameter
    phi = deg2rad(5);
    
    % big radius parameter
    R = 2;
    
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

% this D-curve has sharp edges which cause the arrow in the plot to jump!
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
function g = getTransferFunction(zeros,poles)
    g = @(t) polyval(poly(zeros),t)./polyval(poly(poles),t);
end

%% various stuff
% this function clips the input values to the maximum value allowed inside
% the plot
function z = trunc(x)
    global windowLimits;
    
    z = max(windowLimits.xmin, min(windowLimits.xmax, real(x))) ...
        + i*max(windowLimits.ymin, min(windowLimits.ymax, imag(x)));
end