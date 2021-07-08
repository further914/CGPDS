function result = contains(string, pattern)

    if numel(pattern) == 0
        result = true;
    else
        result = numel(strfind(string, pattern)) ~= 0;
    end

end