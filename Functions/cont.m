function u = cont(phase, x)
    lambda = -x(1,:)';
    phi = 0;
    u = lambda .* (phase~=0) + phi .* (phase==0);
end