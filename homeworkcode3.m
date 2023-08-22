fprintf("Bisections \n");
bsctn(1*10^-5);
bsctn(1*10^-6);
bsctn(1*10^-8);
bsctn(1*10^-10);
fprintf("\n");
fprintf("Fixed Points \n");
fxdPnts(1*10^-5);
fxdPnts(1*10^-6);
fxdPnts(1*10^-8);
fxdPnts(1*10^-10);

function [out] = func(num)
    out = num - cos(num);
end

function [out2] = func2(num)
    out2 = cos(num);
end

function bsctn(tol)
    a = 0;
    b = 1;
    c = 0;
    cntr = 0;
    i = ceil((log(1) - log(tol)) / log(2));
    while(abs(func(b)) > tol)
        c = (a + b) / 2;
        if(func(c) * func(a) < 0)
            b = c;
        else
            a = c;
        end
        cntr = cntr + 1;
    end
    fprintf("Iteration: %f Root: %f \n", i, c);
end

function fxdPnts(tol)
    diff = 10;
    new = 0;
    tempA = [];
    cntr = 1;
    while(abs(diff) > tol)
        new = func2(new);
        tempA(cntr) = new;
        if cntr > 1
            diff = new - tempA(cntr - 1);
        end
        cntr = cntr + 1;
    end
    fprintf("Iteration: %f Root: %f \n", cntr - 1, new);
end
