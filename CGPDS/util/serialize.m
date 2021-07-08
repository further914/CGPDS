function serialize(data, path)

    variable.data = data;
    save(path, '-struct', 'variable');

end