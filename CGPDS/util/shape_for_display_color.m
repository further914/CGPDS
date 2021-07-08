function image = shape_for_display_color(S, H, W)

    L = size(S, 1);
    colors = generate_colors(L);
    image = reshape(S' * colors, [H W 3]);
    
end