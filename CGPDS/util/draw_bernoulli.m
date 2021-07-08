function samples = draw_bernoulli(p)

    samples = double(p > rand(size(p)));

end
