% Makes the plot rectangular (make the smaller axis as large as the bigger one)
function [xlimit,ylimit] = make_limits_quadratic(~, xlimit, ylimit) % ignored argument is 'this'
    
    width = xlimit(2) - xlimit(1);
    height = ylimit(2) - ylimit(1);
    
    if width >= height
        y0 = mean(ylimit);
        ylimit = [(y0-width/2),(y0+width/2)];
    else
        x0 = mean(xlimit);
        xlimit = [(x0-height/2),(x0+height/2)];
    end
    
end
        