estimatedDerivative(.01)
function estimatedDerivative(h)
    hVal = h;
    hVal2 = h/2;
    for i = 0:10
        x = i / 10;
        derivative = 10000*(1/12*hVal)*(getOutput(x-2*hVal)-8*getOutput(x-hVal)+8*getOutput(x+hVal)-getOutput(x+2*hVal));
        derivative2 = 40000*(1/12*hVal2)*(getOutput(x-2*hVal2)-8*getOutput(x-hVal2)+8*getOutput(x+hVal2)-getOutput(x+2*hVal2));
        fprintf("Derivative when x is %i/10 = %f", i, derivative);
        realDer = getOutput(x) + exp(x);
        diff = abs(realDer - derivative);
        diff2 = abs(realDer - derivative2);
        fprintf("    Absolute error when x is %i/10 = %e", i, diff);
        order = log2(diff/diff2);
        fprintf("    Order when x is %i/10 = %f \n", i, order);
    end
end

function output = getOutput(a)
    output = a * exp(a);
end