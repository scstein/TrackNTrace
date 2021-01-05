function d = read_photons(filename, range_seconds)
% Fixed error on too large range. Now the range is constraint to the
% available range. CT, 2020
ms = photonscore.file_read(filename, '/photons/ms');

if nargin == 1
    offset = 0;
    count = ms(end);
else
    ms_range = int32(range_seconds(:) * 1000 + 1);
    ms_range = min(ms_range, length(ms));
    a = min(ms_range);
    b = max(ms_range);
    offset = ms(a);
    count = ms(b) - offset;
end

d.x = photonscore.file_read(filename, '/photons/x', offset, count);
d.y = photonscore.file_read(filename, '/photons/y', offset, count);
d.dt = photonscore.file_read(filename, '/photons/dt', offset, count);
try
    d.t = photonscore.file_read(filename, '/photons/aat', offset, count);
catch e
    warning(e.message)
end
d.ms = ms;
end
