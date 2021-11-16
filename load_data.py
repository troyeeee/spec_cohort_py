def load_geno(path):
    res = []
    with open(path) as f:
        for l in f:
            res.append(l.strip())
    return res