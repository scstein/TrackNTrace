function h = median(fi)
cnt = fi.image(:);
ptr = [0; cumsum(cnt)];
h = zeros(size(fi.image));
for i=1:length(cnt)
    if cnt(i) == 0
        continue
    end
    if mod(cnt(i), 2) == 1
        h(i) = fi.time(ptr(i) + (cnt(i) + 1)/2);
    else
        a = ptr(i) + cnt(i) / 2;
        h(i) = (fi.time(a) + fi.time(a+1)) / 2;
    end
end