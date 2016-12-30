function legacy = anp_check_Matlab_version()
    % release history and numbers:
    % https://en.wikipedia.org/wiki/MATLAB#Release_history
    
	% https://ch.mathworks.com/help/doc-archives.html
	
    % whole release download:
    % https://ch.mathworks.com/downloads/web_downloads/select_release
    
    % only runtime download:
    % https://ch.mathworks.com/products/compiler/mcr/
    
    v = ver;
    idx = find(ismember({v(:).Name},'MATLAB'));
    matlab_version = str2double(v(idx).Version(1:3));

    if matlab_version <= 8.3
        % we want at least HG2 to work with.
        legacy = 'too_old';
        if ~exist('animated_nyquist_plot_warningOff','file')
            fprintf('Warning: You are using a release of Matlab that doesn''t support the current\n')
            fprintf('HG2 graphics engine. Some features such as axle repositioning are deactivated, but there''s still a chance that the script crashes.\n')
            fprintf('Please consider upgrading to R2014b or higher.\n');
            fprintf('This script has not been tested with releases lower than R2014b.\n');
            fprintf('(Suppress this warning by placing an empty file named ''animated_nyquist_plot_warningOff'' in the script''s working folder.)\n');
        end
    elseif matlab_version < 8.6
        legacy = 'pre_R2015b';
    else
        legacy = 'R2015b_or_newer';
    end
    
    if ~exist('animated_nyquist_plot_warningOff','file') && matlab_version >= 8.5 && matlab_version <= 9.0
        fprintf('Note: this script uses release-specific functions of Matlab R2014b\n');
        fprintf('respectively from R2015b on. Testing has been done only for R2014b and R2016b though.\n');
        fprintf('Please report any bugs together with their input!\n');
        fprintf('(Suppress this warning by placing an empty file named ''animated_nyquist_plot_warningOff'' in the script''s working folder.)\n');
    end
end