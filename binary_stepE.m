function output = binary_stepE(x,Edesp)
    if x <= Edesp
        output = 0;
    else
        output = x - Edesp;
    end
end
