function [] = anp_check_Matlab_version()
    %   -------------------------------------------------------------------
    %   Release numbers
    %   ***************
    %   
    %   8.4     R2014b
    %   8.5     R2015a
    %   8.5     R2015aSP1
    %   8.6     R2015b
    %   9.0     R2016a
    %   9.1     R2016b
    %   9.2     R2017a
    %	
    %   -------------------------------------------------------------------
    
    % https://ch.mathworks.com/help/doc-archives.html
	
    % whole release download:
    % https://ch.mathworks.com/downloads/web_downloads/select_release
    
    % only runtime download:
    % https://ch.mathworks.com/products/compiler/mcr/
    
    if any(strcmp(who('global'),'matlab_version'))
        % 'ver' and the following search for the Matlab version string
        % perform badly, so simply reuse the old result.
        global matlab_version;                                  %#ok<*TLEV>
    else
        % 'matlab_version' hasn't been set yet, so continue the version
        % check
        global matlab_version;                                  %#ok<REDEF>
        
        v = ver;
        idx = find(ismember({v(:).Name},'MATLAB'));
        matlab_version = str2double(v(idx).Version(1:3));       %#ok<NASGU>
    end

    if matlab_version < 9.0                                     %#ok<*NODEF>
        % We want at least R2016b to work with in order to use
        % property type validation in classes.
        error('anp:incompatibleRelease', 'You are using a Matlab release prior to R2016a.\nUnfortunately this program needs at least release R2016a to run.\n(Testing has been done on R2016b.)\n');
    elseif matlab_version == 9.0 || matlab_version > 9.1
        if ~exist('anp_warningOff','file')
            warning('anp:unsupportedRelease', 'This program uses release-specific functions of Matlab R2016b.\nTesting has only been done thoroughly R2016b.\nPlease report any bugs together with their input!\n(Suppress this warning by placing an empty file named ''anp_warningOff'' in the program''s working folder.)\n\n');
        end
    end
    
end

