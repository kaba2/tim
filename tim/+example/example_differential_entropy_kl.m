clear all;
close all;

for i = 0 : 2
    for j = 0 : 2
        p = 4 * i + j + 1;
        d = 2^i;
        k = 2^j;
        figure;
        example.draw_differential_entropy_kl('k', k, 'd', d);
    end
end
