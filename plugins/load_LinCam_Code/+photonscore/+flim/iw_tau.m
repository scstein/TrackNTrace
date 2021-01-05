function rgb = iw_tau(a1, a2, a3, a4, a5)
switch nargin
  case 5
    counts = a1;
    counts_range = a2;
    tau = a3;
    tau_range = a4;
    pal = a5;
  case 2
    counts = a1;
    counts_range = auto_range(a1, 100, 10);
    counts_range(1) = 0;
    tau = a2;
    tau_range = auto_range(a2, 5, 5);
    pal = 'preview.png';
  case 3
    counts = a1;
    counts_range = auto_range(a1, 100, 10);
    counts_range(1) = 0;
    tau = a2;
    tau_range = auto_range(a2, 5, 5);
    pal = a3;
  otherwise
    error 'Expecting 2, 3 or 5 arguments';
end

if ischar(pal)
  try
    pal = imread(pal);
  catch
    palpath = [fileparts(mfilename('fullpath')) '/palette/'];
    pal = imread([palpath pal]);
  end
end

counts = double(counts);
counts_range = double(counts_range);
tau = double(tau);
tau_range = double(tau_range);

i1 = x_to_i(counts, counts_range, size(pal, 1));
i2 = x_to_i(tau, tau_range, size(pal, 2));
ix = i1 * size(pal, 2) + i2 + 1;
lin_pal = reshape(pal, [size(pal,1)*size(pal,2) size(pal,3)]);
rgb = reshape(lin_pal(ix, :), [size(counts) size(pal,3)]);
end

function ix = x_to_i(x, xr, n)
ix = floor((x(:) - xr(1)) / (xr(2) - xr(1)) * n);
ix = min(ix, n - 1);
ix = max(ix, 0);
end

function r = auto_range(x, fl, fu)
xx = sort(x(:));
xx = xx(xx > 0);
n1 = floor(length(xx) / fl);
n2 = floor(length(xx) - length(xx) / fu);
n1 = max(1, min(n1, length(xx)));
n2 = max(1, min(n2, length(xx)));
r = [xx(n1) xx(n2)];
end
