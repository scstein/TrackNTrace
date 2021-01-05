function h = decay_from_mask(fi, mask, from, to, bins)
cnt = fi.image(:);
ptr = [0; cumsum(cnt)];
mask = uint32(mask(:));
assert(isequal(length(mask), length(cnt)));
y = uint32(fi.time - fi.time);
for i=1:length(cnt)
    if cnt(i) == 0
        continue
    end
    a = ptr(i) + 1;
    b = ptr(i) + cnt(i);
    y(a:b) = mask(i);
end

m_min = min(mask);
m_max = max(mask);
n_mask = m_max - m_min + 1;

h = photonscore.hist_2d(...
  uint32(fi.time), from, to, bins, ...
  y, m_min, m_max + 1, n_mask);
