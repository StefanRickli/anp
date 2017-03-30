function dbg_out(varargin)
    global debug_text
    if debug_text
        fprintf(varargin{:});
    end
end
