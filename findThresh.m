function [res] = findThresh(omega)
    thresh = deg2rad(3);
    res = zeros(3,1);
    for j = 1:length(omega)
        if omega(1,j) < thresh
            res(1) = j;
        elseif omega(2,j) < thresh
            res(2) = j;
        elseif omega(3,j) < thresh
            res(3) = j;
        end
        res = max(res);
    end
end