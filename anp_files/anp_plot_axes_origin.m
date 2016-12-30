% moves the current ax(l)es into the origin of the plot
function anp_plot_axes_origin(legacy)
    switch(legacy)
        case 'R2015b_or_newer'
            ax = gca;
            ax.XAxisLocation = 'origin';
            ax.YAxisLocation = 'origin';
        case 'pre_R2015b'
            % http://undocumentedmatlab.com/blog/customizing-axes-part-2
            ax = gca;
            ax.XBaseline.Color = 'k';
            ax.YBaseline.Color = 'k';
            
            ax.XBaseline.LineWidth = 1.5;
            ax.YBaseline.LineWidth = 1.5;
            
            ax.YBaseline.Visible = 'on';
            ax.XBaseline.Visible = 'on';
            
            ax.XRuler.Axle.Visible = 'off';
            ax.YRuler.Axle.Visible = 'off';

            ax.LineWidth = 1;

            % 2nd possibility for workaround: http://matlabvn.com/matlab-plotting-tricks-move-x-axis-center/
            % other:
            % https://ch.mathworks.com/matlabcentral/answers/147087-how-to-change-the-axes-position-in-matlab
            % http://ch.mathworks.com/matlabcentral/fileexchange/22956-axescenter
        otherwise
            % Too old to produce compatible code. Just skip the axle
            % repositioning.
            return;
    end
end