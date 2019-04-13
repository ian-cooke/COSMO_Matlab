function [m] = bangControl(bp)
    global m_max;
     m = -m_max*sign(bp);
end

