function output = binary_stepI(x,Idesp)
    if x <= Idesp
        output = 0;
    else
        output = x - Idesp;
    end
end
