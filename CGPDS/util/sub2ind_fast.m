function ind = sub2ind_fast(dimensions, i, j)

    ind = (j-1) * dimensions(1) + i;

end