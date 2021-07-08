function data = deserialize(path)

    tmp = load(path);
    data = tmp.data;

end