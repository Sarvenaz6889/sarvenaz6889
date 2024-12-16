function sig = glrtSignal(params, timeVec)
    a1 = params(1);
    a2 = params(2);
    a3 = params(3);
    sig = sin(2 * pi * (a1 * timeVec + a2 * timeVec.^2 + a3 * timeVec.^3));
end
