function y = repmat_fast(x, repeat_y)

    y = x(ones(1, repeat_y), :);
    
end