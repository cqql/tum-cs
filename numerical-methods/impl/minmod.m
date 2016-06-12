function res = minmod(a, b)
    if a*b <0
        res = 0;
    elseif abs(a)<abs(b)
        res = a;
    else 
        res = b;
    end
end

