function psum = partials(y)
[a,b] = size(y);
psum = zeros(a,b);
psum(1) = y(1);
for i = 2:a
    psum(i) = psum(i - 1) + y(i);
end