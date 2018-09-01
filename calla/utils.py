from calla import abacus

def table(calculators:[abacus]):
    if calculators == None or len(calculators)<1:
        return None
    types = []
    calc0 = calculators[0]
    attrs = list(calc0.__inputs__.keys())+list(calc0.__deriveds__.keys())
    for i in range(1,len(calculators)):
        calc = calculators[i]
        t = type(calc)
        if not t in types:
            types.append(t)
            attrs_t = list(t.__inputs__.keys())+list(t.__deriveds__.keys())
            attrs = [attr for attr in attrs if attr in attrs_t]
    # build table
    tbl = []
    tbl.append([calc0.symbol(attr) for attr in attrs])
    for calc in calculators:
        paras = calc.parameters()
        tbl.append([paras[attr] for attr in attrs])
    return tbl

