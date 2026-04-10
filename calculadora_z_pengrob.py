
import matplotlib.pyplot as plt 
import numpy as np 

print(f"---------------------------------")
print(f"Cálculo de Z usando Peng-Robinson")
print(f"---------------------------------")
print(f"------------Dados-----------------") 

# Dados:
esp = input("Nome da Substância: ")
Tc = float(input("Temperatura Crítica [K] : "))
Pc = float(input("Pressão Crítica [bar]: "))*1e5 
w = float(input("Fator acêntrico da espécie:"))
T = float(input("Temperatura de operação [K]:"))
P = float(input("Pressão de operação [bar]:"))*1e5
R = 8.314 
oma = 0.45724
omb = 0.07780

# Variáveis:
Tr = (T/Tc)
kappa = (0.37464 + 1.54226*w - 0.26992*(w**2))
alfa = (1 + kappa*(1-(Tr)**0.5))**2

# Definindo os parâmetros:
def a(Tc,Pc): 
    res = oma*((R**2)*(Tc**2)/Pc)
    return res*alfa

def b(Tc, Pc):
    res2 = omb*(R*Tc/Pc)
    return res2
 
A = a(Tc,Pc)*P/((R**2)*(T**2)) 
B = b(Tc,Pc)*P/(R*T)


def f(z):
    return (z**3)-((1 - B)*z**2)+((A - 2*B - 3*B**2)*z)-(A*B - B**2 - B**3) 

def der(f, z, h=1e-7):
    return (f(z+h) - f(z-h))/(2*h)


# Cálculo da iteração: 
itmax = 1000
tol = 1e-6 
raizes = []

for z0 in [B, 0.1, 1.0]:
    f0 = f(z0) 
    df0 = der(f, z0)
    passo = -f0/df0
    IT = 0 

    while abs(f0) > tol and IT < itmax:
        z1 = z0 + passo
        f1 = f(z1)
        
        if abs(f1) < abs(f0):
            z0 = z1 
            f0 = f1
            df0 = der(f, z0)
            passo = -f0/df0
        else:
            passo = passo/2 
        
        IT +=1

    if z0 > B:
        raizes.append(round(z0, 6))

raizes = sorted(set(raizes))

# Resultados:  
print(f"==================================")
print(f"-------------RESULTADO------------")
print(f"==================================")
print(f"Espécie Química: {esp}")
print(f"Solução:    Z = {z0:.6f}")
print(f"Iterações:  {IT}")
print(f"Resíduo:    {abs(f0):.2e}") 

if len(raizes) == 1:
    print(f"Z = {raizes[0]:.6f} (fase única)") 
   

elif len(raizes) == 2:
    print(f"Z líquido = {raizes[0]:.6f}")
    print(f"Z vapor   = {raizes[1]:.6f}")
    print(f"Região bifásica — use fugacidade para determinar fase estável")

elif len(raizes) == 3:
    print(f"Z líquido      = {raizes[0]:.6f}")
    print(f"Z intermediário = {raizes[1]:.6f} (instável, descartado)")
    print(f"Z vapor        = {raizes[2]:.6f}") 

# Plotagem do gráfico:
z = np.linspace(0.01, 1.50, 500)
plt.figure(figsize=(15,10))
plt.plot(z, f(z)) 
plt.axhline(0, color='black', linewidth=1.5)
plt.title(f"Raízes Cúbicas de Z - {esp} | Temp. Operação: {T} [K] | Pressão Operação: {P/1e5} [bar]")
plt.xlabel("Z")
plt.ylabel("f(Z)")

if len(raizes) == 1:
    plt.plot(raizes[0], f(raizes[0]), 'o', color='red', markersize=8, label=f'Z - Líquido = {raizes[0]:.4f}')  

elif len(raizes) == 2:
    plt.plot(raizes[1], f(raizes[1]), 'o', color='blue', markersize=8, label=f'Z - Intermediário = {raizes[1]:.4f}')
    plt.plot(raizes[2], f(raizes[2]), 'o', color='purple', markersize=8, label=f'Z - Vapor = {raizes[2]:.4f}')

elif len(raizes) == 3:
    plt.plot(raizes[0], f(raizes[0]), 'o', color='red', markersize=8, label=f'Z - Líquido = {raizes[0]:.4f}')
    plt.plot(raizes[1], f(raizes[1]), 'o', color='blue', markersize=8, label=f'Z - Intermediário = {raizes[1]:.4f}')
    plt.plot(raizes[2], f(raizes[2]), 'o', color='purple', markersize=8, label=f'Z - Vapor = {raizes[2]:.4f}')

plt.legend(loc='best', shadow=True, fontsize='medium')
plt.grid(True, linestyle='--', alpha=0.7)
plt.show()