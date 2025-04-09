syms x
g=2/x;
for i = 1:7
    N(i,1) = double((subs(g, c+h) - subs(g, c-h)) / (2*h));
    h = h / 2;
end

for j = 2:7
    for i = 1:7-j+1
        N(i,j) = (4^(j-1) * N(i+1,j-1) - N(i,j-1)) / (4^(j-1) - 1);
    end
end

N(1,7)
