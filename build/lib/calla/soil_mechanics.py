def soil_saturation_density(Gs, e, water_density=10):
    return (Gs+e)/(1+e)*water_density

def test():
    Gs = 2.7
    e = 0.744
    rho = soil_saturation_density(Gs, e)
    print(rho)

if __name__ == '__main__':
    test()
