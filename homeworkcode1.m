fprintf("Sum of squares in matrix: %f\n", sqSum([12.2 1.2 2.4; 2 3 4; 2.5 6.2 3.4]));
fprintf("20.25 to bin: %s\n", decToBin(20.25));
fprintf("21.45 to bin: %s\n", decToBin(21.45));

function [runningSum] = sqSum(test)
    runningSum = 0;
    matSize = size(test);
    for r = 1:matSize(1)
        for c = 1:matSize(2)
            runningSum = runningSum + test(r,c)^2;
        end
    end
end

function [binConv] = decToBin(decNum)
    temp = "";
    binConv = "";
    if (decNum < 0)
        decNum = -1 * decNum;
        binConv = "-";
    end
    heldQ = floor(decNum);
    while heldQ ~= 0
        rem = mod(heldQ, 2);
        heldQ = floor(heldQ/2);
        temp = strcat(temp, string(rem));
    end
    binInt = reverse(temp);
    cntr = 0;
    lastBits = "";
    fraction = abs(decNum - floor(decNum));
    while cntr ~= 10 && fraction ~= 0
        fraction = fraction * 2;
        lastBits = lastBits + string(floor(fraction));
        fraction = abs(fraction - floor(fraction));
        cntr = cntr + 1;
    end
        binConv = binConv + binInt + "." + lastBits;
end