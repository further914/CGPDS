function colors = generate_colors(L)

    colors = jet;
    colors = colors(1:floor(63/(L-1)):64, :);

    if L == 1
        colors = [1 1 1];
    elseif L == 2
        colors = [1 1 1; 0 0 0];
    elseif L == 3
        colors = [0 0 0; 1 0.8 0;  0.75 0.3 0;];
    elseif L == 7
        colors = [       0         0         0
                         0         0    1.0000
                         0    1.0000         0
                         0    1.0000    1.0000
                    0.5020    0.5020    0.5020
                    1.0000         0         0
                    1.0000    1.0000         0];
    elseif L == 8
        colors = [      0         0    0.5882
                         0    0.1765    1.0000
                         0    0.7647    1.0000
                    0.3529    1.0000    0.6471
                    0.6471         0         0
                    0.8824         0         0
                    0.9412    1.0000    0.0588
                    1.0000    0.4706         0];
    end

end