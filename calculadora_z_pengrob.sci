clc; clear 

// Cálculo de Z usando Peng-Robinson

// @Author: @Nogfe4

disp("---------------------------------")
disp("Cálculo de Z usando Peng-Robinson")
disp("---------------------------------")
disp("------------Dados-----------------")

// Dados:
esp = input("Nome da Substância: ", "string")
Tc  = input("Temperatura Crítica [K] : ")
Pc  = input("Pressão Crítica [bar]: ") * 1e5
w   = input("Fator acêntrico da espécie: ")
T   = input("Temperatura de operação [K]: ")
P   = input("Pressão de operação [bar]: ") * 1e5

R   = 8.314
oma = 0.45724
omb = 0.07780

// Variáveis:
Tr    = T / Tc
kappa = 0.37464 + 1.54226*w - 0.26992*(w^2)
alfa  = (1 + kappa*(1 - Tr^0.5))^2

// Definindo os parâmetros:
function res = a(Tc, Pc)
    res = oma * ((R^2) * (Tc^2) / Pc) * alfa
endfunction

function res2 = b(Tc, Pc)
    res2 = omb * (R * Tc / Pc)
endfunction

A = a(Tc, Pc) * P / ((R^2) * (T^2))
B = b(Tc, Pc) * P / (R * T)

function val = f(z)
    val = (z^3) - ((1-B)*z^2) + ((A - 2*B - 3*B^2)*z) - (A*B - B^2 - B^3)
endfunction

function val = der(z)
    h   = 1e-7
    val = (f(z+h) - f(z-h)) / (2*h)
endfunction

// Cálculo da iteração:
itmax  = 1000
tol    = 1e-6
raizes = []

chutes = [B, 0.1, 1.0]

for i = 1:3
    z0    = chutes(i)
    f0    = f(z0)
    passo = -f0 / der(z0)
    IT    = 0

    while abs(f0) > tol & IT < itmax
        z1 = z0 + passo
        f1 = f(z1)

        if abs(f1) < abs(f0) then
            z0    = z1
            f0    = f1
            passo = -f0 / der(z0)
        else
            passo = passo / 2
        end

        IT = IT + 1
    end

    if z0 > B then
        raizes = [raizes, round(z0 * 1e6) / 1e6]
    end
end

raizes = unique(raizes)
raizes = gsort(raizes, 'g', 'i')
nr     = length(raizes)

// Resultados:
disp("==================================")
disp("-------------RESULTADO------------")
disp("==================================")
disp("Espécie Química: " + esp)
disp("Solução:    Z = " + string(z0))
disp("Iterações:  " + string(IT))
disp("Resíduo:    " + string(abs(f0)))

if nr == 1 then
    disp("Z = " + string(raizes(1)) + " (fase única)")

elseif nr == 2 then
    disp("Z líquido = " + string(raizes(1)))
    disp("Z vapor   = " + string(raizes(2)))
    disp("Região bifásica — use fugacidade para determinar fase estável")

elseif nr == 3 then
    disp("Z líquido       = " + string(raizes(1)))
    disp("Z intermediário = " + string(raizes(2)) + " (instável, descartado)")
    disp("Z vapor         = " + string(raizes(3)))
end

// Plotagem do gráfico:
z_arr = linspace(0.01, 1.50, 500)
f_arr = zeros(1, 500)
for i = 1:500
    f_arr(i) = f(z_arr(i))
end

scf()
plot(z_arr, f_arr)
xgrid()
plot([0.01, 1.50], [0, 0], 'k-')
xtitle("Raízes Cúbicas de Z - " + esp + " | Temp. Operação: " + string(T) + " [K] | Pressão Operação: " + string(P/1e5) + " [bar]", "Z", "f(Z)")

cores = ['red', 'blue', 'purple']

if nr == 1 then
    plot(raizes(1), f(raizes(1)), 'ro', 'MarkerSize', 8)
    legend("Z - Líquido = " + string(raizes(1)))

elseif nr == 2 then
    plot(raizes(1), f(raizes(1)), 'bo', 'MarkerSize', 8)
    plot(raizes(2), f(raizes(2)), 'mo', 'MarkerSize', 8)
    legend(["Z - Intermediário = " + string(raizes(1)), "Z - Vapor = " + string(raizes(2))])

elseif nr == 3 then
    plot(raizes(1), f(raizes(1)), 'ro', 'MarkerSize', 8)
    plot(raizes(2), f(raizes(2)), 'bo', 'MarkerSize', 8)
    plot(raizes(3), f(raizes(3)), 'mo', 'MarkerSize', 8)
    legend(["Z - Líquido = " + string(raizes(1)), "Z - Intermediário = " + string(raizes(2)), "Z - Vapor = " + string(raizes(3))])
end
