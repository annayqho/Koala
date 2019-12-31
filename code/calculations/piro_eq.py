E = 10
K = 0.1

# First, solve for R_e
Re = (1.9/(K/0.34)**(-1.4))**4 * 1E13
print(Re/1E11)

# Then solve for the core mass
Mc = 3 * ((Re/1.4E14) * (1/(K/0.34)**0.74) * (1/E**(-0.87)))**(1/0.61)
print(Mc)

# Then solve for the envelope mass
Me = 0.027 * (E)**(0.43) * (K/0.34)**(-0.87) * (Mc/3)**(-0.3)
print(Me)
