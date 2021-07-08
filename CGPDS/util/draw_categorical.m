function draw = draw_categorical(p)

    % get the data dimensions
    [L, D] = size(p);
    
    % draw
    cumulative_p = cumsum(p, 1);
    hit_p = cumulative_p > repmat_fast(rand(1, D), L);
    draw_categories = (L+1)-sum(hit_p,1);
    
    % sometimes you end up with dimensions which are impossible, i.e. every
    % category is impossible. in these cases, just choose on category at
    % random.
    impossibles = find(draw_categories==L+1);
    draw_categories(impossibles) = floor(rand(1,numel(impossibles))*L)+1;
    
    % create the output matrix
    draw_indices = sub2ind_fast([L,D],draw_categories,(1:D));
    draw = zeros(L,D);
    draw(draw_indices) = 1;
    
end