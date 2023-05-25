function y = soft(x,T)
if length(T)==1
T = T + eps;
y = max(abs(x) - T, 0);
y = y./(y+T) .* x;
else 
    T=repmat(T',size(x,1),1);
    T = T + eps;
    y = max(abs(x) - T, 0);
    y = y./(y+T) .* x;
end

