function [y, outs] = eval_ann(w, layers, x)
w = w(:);
layers = [size(x, 1); layers(:); 1];
n = length(layers);
wp = 1;
y = x;
outs = cell(1, n-1);
for i=1:n-1
    ins = layers(i);
    ots = layers(i + 1);
    B = w(wp:wp+ots-1);
    wp = wp + ots;
    W = reshape(w(wp:wp+ins*ots-1), ots, ins);
    wp = wp + ots * ins;
    y = W*y + repmat(B, 1, size(y, 2));
    if i ~= n-1
        outs{i} = y;
        y = tanh(y);
    end
end

