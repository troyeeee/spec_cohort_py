def comp(in_str):
    res = ""
    for s in in_str:
        if s=="0":
            res=res+"1"
        elif s =="1":
            res=res+"0"
        else:
            res+=s
    return res