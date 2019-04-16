function Gamma = calc_gamma(Bb_prev, Bb, t)
%CALC_GAMMA calcualate matrix for lovera track control

global tstep;

S_prev = tilde(Bb_prev./norm(Bb_prev));

S = tilde(Bb./norm(Bb));


prod_prev = S_prev*S_prev';
prod = S*S';

Gamma = zeros(3,3);
for i = 1:3
    for j = 1:3
        Gamma(i,j) = trapz([t, t+tstep], [prod_prev(i,j), prod(i,j)]);
    end
end

end