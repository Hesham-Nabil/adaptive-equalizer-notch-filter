function [yn, xn1] = adaptive_filter_OH(xn, an, yn, u, r, k)
    if k == 1
        en = xn(k);
        yn(k) = en;
        xn1 = 0;
    elseif k == 2
        en = xn(k) + an(k) * xn(k-1);
        yn(k) = en - r*an(k)*yn(k-1);
        xn1 = xn(k-1);
    else
        en = xn(k) + an(k)*xn(k-1) + xn(k-2);
        yn(k) = en - r*an(k)*yn(k-1) - r^2*yn(k-2);
        xn1 = xn(k-1);
    end
end