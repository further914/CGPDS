function p = softmax(q)

    r = exp(q);
    p = r ./ repmat_fast(sum(r,1), size(r,1));
    
    if any(any(isnan(p)))
        
        % revert to the proper way
        rows = size(q, 1);
        r = exp(q - repmat_fast(max(q,[],1), rows) + log(realmax)/2);
        p = r ./ repmat_fast(sum(r,1), rows);
        
    end
    
end