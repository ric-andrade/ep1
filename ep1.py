import math
import matplotlib.pyplot as plt
import numpy as np
# import sys
np.set_printoptions(precision=4)  # Limitar precisão dos dados. Obs.: Só afeta os prints relacionados ao Numpy.


def main():
    """Função principal responsável por chamar as demais funções contidas nesse programa,
    executando, assim, o algoritmo."""

    global n    # Definição da variável global "n" que será utilizada em diversas funções neste programa.
    global gr   # Definição da variável global "gr" que será utilizada como parâmetro para plotar ou não os gráficos.
    met = input("Qual método deseja utilizar? ('d', 'e' ou 'k') ")
    # t = float(input("Digite o valor de tempo T: "))
    t = 1
    n = int(input("Digite o valor de N: "))
    if met != "d":
        lda = n             # Devido à condição (Delta t = Delta x)
        m = round(n * t)
    else:
        lda = float(input("Digite o valor de lambda: "))
        m = (t*(n**2))/lda

        if m != int(m):          # Arredondamento do parâmetro "m" segue convenção explicada em aula.
            decimal = m - int(m)
            arredondado = int(m)

            if decimal >= 0.5:
                arredondado += 1
                m = arredondado

            else:
                m = arredondado

        else:
            m = int(m)

    teste = input("Qual teste deseja realizar? ('a', 'b' ou 'c') ")
    gr = input("Deseja visualizar os gráficos? (y: yes ; n: no) ")
    v = crie_vetor(n + 1)
    v_contorno = contorno(v, teste)       # Condição inicial (t = 0)
    v_exato_contorno = v.copy()
    v_exato_contorno[0] = g1(t, teste)
    v_exato_contorno[n] = g2(t, teste)
    if met != "d":
        v_aprox = u(v_contorno, m, t, lda, met, teste)
    else:
        v_aprox = u(v_contorno, m, t, lda, met, teste)

    print("O vetor aproximado no instante tm para cada posição xi é: \n")
    print(np.array(v_aprox))

    if teste != "c":
        v_exato = vetor_exato(v_exato_contorno, m, t, m, teste)
        print("O vetor de resultados exatos no instante tm para cada posição xi é: \n")
        print(np.array(v_exato))
        print("O vetor dos erros, em módulo, no instante tm para cada posição xi é: \n")
        vet_erro = []
        for i in range(n + 1):
            vet_erro.append(abs(v_aprox[i] - v_exato[i]))
        print(np.array(vet_erro))
        print("O erro máximo no instante tm é: ")
        print(max(vet_erro))


def plot(v_aprox, v_exato, teste):
    """ Função "plot". Recebe três argumentos, sendo responsável por plotar os gráficos deste programa.

    Arguments:
        v_aprox {list} -- vetor de temperaturas aproximadas
        v_exato {list} -- vetor de temperaturas exatas esperadas
        teste {string} -- variável que contém o tipo de teste que o usuário deseja executar
    """    

    if teste != "c":
        position = []               # Plot do gráfico para os itens A e B.
        for i in range(n + 1):
            position.append((1 / n) * i)

        plt.scatter(position, v_aprox, marker='*', c='r')
        plt.scatter(position, v_exato, marker='o', c='b')
        plt.scatter(position, np.abs(np.array(v_aprox) - np.array(v_exato)), marker='^', c='k')

        plt.xlabel('position')
        plt.ylabel('temperature')
        plt.title('Temperature X Position')
        plt.legend(['aproximate solution', 'exact solution', 'error'])
        plt.grid(True)
        plt.show()

    else:
        position = []               # Plot do gráfico para o item C.
        for i in range(n + 1):
            position.append((1 / n) * i)

        plt.scatter(position, v_aprox, marker='*', c='r')

        plt.xlabel('position')
        plt.ylabel('temperature')
        plt.title('Temperature X Position')
        plt.legend(['aproximate solution'])
        plt.grid(True)
        plt.show()

    if teste != "c" and gr == "y":
        vet_erro = []
        position = []                   # plot do gráfico dos erros.
        for i in range(n + 1):
            vet_erro.append(abs(v_aprox[i] - v_exato[i]))
            position.append((1 / n) * i)

        plt.scatter(position, vet_erro, marker='^', c='k')

        plt.xlabel('position')
        plt.ylabel('error')
        plt.title('Error X Position')
        plt.legend(['error'])
        plt.grid(True)
        plt.show()


def decomp(vet1, vet2):
    """ Função "decomp". Recebe dois argumentos e é responsável por decompor a matriz A em outras três matrizes (L, D e L^t) que podem ser representadas por dois vetores "vet_d" e "vet_l". A matriz A deve ser tridiagonal simétrica contendo valores nulos em posições diferentes da diagonal principal e das duas subdiagonais. A decomposição respeita a seguinte relação: A = LDL^t.

    Arguments:
        vet1 {list} -- vetor que representa a diagonal principal da matriz A
        vet2 {list} -- vetor que representa as duas subdiagonais da matriz A

    Returns:
        list -- lista composta pelos vetores "vet_d" e "vet_l" que representam a matriz A decomposta
    """    

    vet_d = crie_vetor(len(vet1))
    vet_l = crie_vetor(len(vet2))

    vet_d[0] = vet1[0]

    h = 0
    while h < n - 2:
        vet_l[h] = vet2[h] / vet_d[h]
        vet_d[h + 1] = vet1[h + 1] - (vet_l[h] ** 2) * vet_d[h]

        h += 1

    return [vet_d, vet_l]


def solve_sys(unew, vet_l, vet_d, vet_dir):
    """ Função "solve_sys". Recebe quatro argumentos que caracterizam um sistema do tipo LDL^t x = b, e retorna o vetor solução "x" deste sistema.

    Arguments:
        unew {list} -- vetor incógnita do sistema linear "x"
        vet_l {list} -- vetor que representa a matriz L
        vet_d {list} -- vetor que representa a matriz D
        vet_dir {list} -- vetor coluna que representa o lado direito da equação matricial "b"

    Returns:
        list -- vetor incógnita do sistema linear "x"
    """    

    y = crie_vetor(n - 1)  # Resolução do sistema Ly = b

    y[0] = vet_dir[1]

    for i in range(1, n - 1):
        y[i] = vet_dir[i + 1] - vet_l[i - 1] * y[i - 1]

    z = crie_vetor(n - 1)  # Resolução do sistema Dz = y

    for o in range(n - 1):
        z[o] = y[o] / vet_d[o]

    unew[n - 1] = z[n - 2]

    for q in range(n - 2, 0, -1):  # Resolução do sistema L^t x = z

        unew[q] = z[q - 1] - vet_l[q - 1] * unew[q + 1]

    return unew


def sol_exata(t, x, teste):
    """ Função "sol_exata". Recebe três argumentos e retorna o valor da solução exata para determinado teste.

    Arguments:
        t {float} -- variável temporal
        x {float} -- variável espacial
        teste {string} -- variável que contém o tipo de teste que o usuário deseja executar

    Returns:
        float -- valor da solução exata correspondente ao teste que o usuário deseja executar 
    """    

    if teste == "a":
        # return 10*t*x**2 * (x - 1)
        return (1 + math.sin(10*t))*(x**2)*((1-x)**2)

    elif teste == "b":
        return math.exp(t-x) * math.cos(5*t*x)


def vetor_exato(vetor, m, t, k, teste):
    """ Função "vetor_exato". Recebe cinco argumentos e retorna o vetor exato da solução da equação do calor para um instante tk.

    Arguments:
        vetor {list} -- vetor que será preenchido
        m {float} -- parâmetro do problema
        t {float} -- variável temporal
        k {int} -- parâmetro do problema
        teste {string} -- variável que contém o tipo de teste que o usuário deseja executar

    Returns:
        list -- vetor exato para determinado instante de tempo tk
    """    

    delta_x = 1 / n
    delta_t = t / m
    i = 1
    while i < n:
        vetor[i] = sol_exata(k * delta_t, i * delta_x, teste)
        i += 1

    return vetor


def fonte(x, t, teste):
    """ Função "fonte". Recebe três argumentos e retorna o valor da função fonte correspondente ao teste que o usuário deseja executar.

    Arguments:
        x {float} -- variável espacial
        t {float} -- variável temporal
        teste {string} -- variável que contém o tipo de teste que o usuário deseja executar

    Returns:
        float -- valor da função fonte correspondente ao teste que o usuário deseja executar 
    """    

    p = 0.25
    delta_x = 1 / n

    if teste == "a":
        # return 10*((x**2)*(x - 1)) - 60*x*t + 20*t
        return 10*math.cos(10*t)*(x**2)*(1-x)**2 - (1+math.sin(10*t))*(12*(x**2)-12*x + 2)

    elif teste == "c":
        if p + (delta_x/2) >= x >= p - (delta_x/2):
            gh = 1/delta_x
        else:
            gh = 0
        return (10000*(1-2*(t**2)))*gh

    else:
        return math.exp(t-x) * (25*t**2 * math.cos(5*t*x) - 10*t * math.sin(5*t*x) - 5*x * math.sin(5*t*x))


def u(vetor, m, t, lda, met, teste):
    """ Função "u". Recebe seis argumentos e é responsável por executar algum método determinado pelo usuário, retornando o vetor aproximado dos valores das temperaturas no instante tM, isto é, no útlimo instante, referente à equação do calor. 

    Arguments:
        vetor {list} -- vetor que contém a condição de contorno u0
        m {float} -- parâmetro do problema
        t {float} -- intervalo de tempo desejado "T"
        lda {float} -- parâmetro do problema
        met {string} -- variável que determina qual o método que o usuário deseja executar
        teste {string} -- variável que contém o tipo de teste que o usuário deseja executar

    Returns:
        list -- vetor de resultados aproximados para o instante tM
    """    

    uold = vetor.copy()          # Recebe o parâmetro "vetor" que contém a condição de contorno u0 (primeira linha)
    unew = crie_vetor(len(vetor))
    delta_x = 1 / n
    delta_t = t / m

    if met == "d":               # Método da discretização (explícito)
        q = 1
        for k in range(m):
            i = 1
            while i < n:
                unew[i] = uold[i] + delta_t * ((uold[i - 1] - 2*uold[i] + uold[i + 1])/delta_x**2 + fonte(i*delta_x, k*delta_t, teste))
                i += 1
            unew[0] = g1((k + 1) * delta_t, teste)
            unew[n] = g2((k + 1) * delta_t, teste)
            uold = np.array(unew.copy())

            if gr == "y":

                if k == q * round((m / 10)):
                    v = crie_vetor(n + 1)
                    v_exato_contorno = v.copy()
                    v_exato_contorno[0] = g1(k * delta_t, teste)
                    v_exato_contorno[n] = g2(k * delta_t, teste)
                    v_exato = vetor_exato(v_exato_contorno, m, t, k, teste)

                    plot(unew, v_exato, teste)
                    q += 1

                elif k == (m - 1):
                    v = crie_vetor(n + 1)
                    v_exato_contorno = v.copy()
                    v_exato_contorno[0] = g1((k+1) * delta_t, teste)
                    v_exato_contorno[n] = g2((k+1) * delta_t, teste)
                    v_exato = vetor_exato(v_exato_contorno, m, t, k+1, teste)

                    plot(unew, v_exato, teste)

    elif met == "e":                    # Método de Euler
        vet1 = []
        for i in range(2, len(vetor)):
            vet1.append(1 + 2 * lda)

        vet2 = []
        for j in range(2, len(vetor) - 1):
            vet2.append(-lda)

        vet_decomp = decomp(vet1, vet2)

        vet_d = vet_decomp[0]
        vet_l = vet_decomp[1]

        vet_dir = crie_vetor(n + 1)     # Início do método de Euler implícito

        r = 1
        for k in range(m):

            for d in range(1, n):       # Preenchimento da matriz coluna b para o método de Euler

                if d == 1:
                    vet_dir[d] = uold[d] + (delta_t * fonte(d*delta_x, (k + 1) * delta_t, teste)) + \
                                 lda * g1((k + 1) * delta_t, teste)
                elif d == n - 1:
                    vet_dir[d] = uold[d] + (delta_t * fonte(d*delta_x, (k + 1) * delta_t, teste)) + \
                                 lda * g2((k + 1) * delta_t, teste)
                else:
                    vet_dir[d] = uold[d] + (delta_t*fonte(d*delta_x, (k + 1)*delta_t, teste))

            unew = solve_sys(unew, vet_l, vet_d, vet_dir)

            unew[0] = g1((k+1)*delta_t, teste)  # Preenchimento da primeira e última posição.
            unew[n] = g2((k+1)*delta_t, teste)
            uold = np.array(unew.copy())

            if gr == "y":

                if k == r * round((m / 10)):
                    v = crie_vetor(n + 1)
                    v_exato_contorno = v.copy()
                    v_exato_contorno[0] = g1(k * delta_t, teste)
                    v_exato_contorno[n] = g2(k * delta_t, teste)
                    v_exato = vetor_exato(v_exato_contorno, m, t, k, teste)

                    plot(unew, v_exato, teste)
                    r += 1

                elif k == (m - 1):
                    v = crie_vetor(n + 1)
                    v_exato_contorno = v.copy()
                    v_exato_contorno[0] = g1((k+1) * delta_t, teste)
                    v_exato_contorno[n] = g2((k+1) * delta_t, teste)
                    v_exato = vetor_exato(v_exato_contorno, m, t, k+1, teste)

                    plot(unew, v_exato, teste)

    else:                               # Método de Crank-Nicolson
        vet1 = []
        for i in range(2, len(vetor)):
            vet1.append(1 + lda)

        vet2 = []
        for j in range(2, len(vetor) - 1):
            vet2.append(-lda/2)

        vet_decomp = decomp(vet1, vet2)

        vet_d = vet_decomp[0]
        vet_l = vet_decomp[1]

        vet_dir = crie_vetor(n + 1)  # Início do método de Crank-Nicolson implícito

        r = 1
        for k in range(m):

            for d in range(1, n):  # Preenchimento da matriz coluna b para o método de Crank-Nicolson

                if d == 1:
                    vet_dir[d] = uold[d] + (lda/2)*(uold[d-1]-2*uold[d]+uold[d+1]) + (delta_t/2) * \
                                 (fonte(d * delta_x, (k + 1) * delta_t, teste) + fonte(d * delta_x, k * delta_t, teste)
                                  ) + (lda/2) * g1((k + 1) * delta_t, teste)
                elif d == n - 1:
                    vet_dir[d] = uold[d] + (lda/2)*(uold[d-1]-2*uold[d]+uold[d+1]) + (delta_t/2) * \
                                 (fonte(d * delta_x, (k + 1) * delta_t, teste) + fonte(d * delta_x, k * delta_t, teste)
                                  ) + (lda/2) * g2((k + 1) * delta_t, teste)
                else:
                    vet_dir[d] = uold[d] + (lda/2)*(uold[d-1]-2*uold[d]+uold[d+1]) + (delta_t/2) * \
                                 (fonte(d * delta_x, (k + 1) * delta_t, teste) + fonte(d * delta_x, k * delta_t, teste))

            unew = solve_sys(unew, vet_l, vet_d, vet_dir)

            unew[0] = g1((k + 1) * delta_t, teste)  # Preenchimento da primeira e última posição.
            unew[n] = g2((k + 1) * delta_t, teste)
            uold = np.array(unew.copy())

            if gr == "y":

                if k == r * round((m / 10)):
                    v = crie_vetor(n + 1)
                    v_exato_contorno = v.copy()
                    v_exato_contorno[0] = g1(k * delta_t, teste)
                    v_exato_contorno[n] = g2(k * delta_t, teste)
                    v_exato = vetor_exato(v_exato_contorno, m, t, k, teste)

                    plot(unew, v_exato, teste)
                    r += 1

                elif k == (m - 1):
                    v = crie_vetor(n + 1)
                    v_exato_contorno = v.copy()
                    v_exato_contorno[0] = g1((k + 1) * delta_t, teste)
                    v_exato_contorno[n] = g2((k + 1) * delta_t, teste)
                    v_exato = vetor_exato(v_exato_contorno, m, t, k + 1, teste)

                    plot(unew, v_exato, teste)

    return unew


def u0(x, teste):
    """ Função "u0". Recebe dois parâmetros e retorna o valor da temperatura correspondente ao instante t=0.

    Arguments:
        x {float} -- variável espacial
        teste {string} -- variável que contém o tipo de teste que o usuário deseja executar

    Returns:
        float -- valor da temperatura correspondente à condição inicial temporal (t=0)
    """    

    if teste == "a":
        # return 0 * x
        return (x**2)*(1-x)**2

    elif teste == "c":
        return 0 * x

    else:
        return math.exp(-x)


def g1(t, teste):
    """ Função "g1". Recebe dois parâmetros e retorna o valor da temperatura correspondente à posição x=0.

    Arguments:
        t {float} -- variável temporal
        teste {string} -- variável que contém o tipo de teste que o usuário deseja executar

    Returns:
        float -- valor da temperatura correspondente à condição inicial espacial (x=0)
    """    

    if teste == "a":
        # return 0 * t
        return 0 * t

    elif teste == "c":
        return 0 * t

    else:
        return math.exp(t)


def g2(t, teste):
    """ Função "g2". Recebe dois parâmetros e retorna o valor da temperatura correspondente à posição x=L, onde L é a posição final barra (L=1).

    Arguments:
        t {float} -- variável temporal
        teste {string} -- variável que contém o tipo de teste que o usuário deseja executar

    Returns:
        float -- valor da temperatura correspondente à condição final espacial (x=1)
    """    

    if teste == "a":
        # return 0 * t
        return 0 * t

    elif teste == "c":
        return 0 * t

    else:
        return math.exp(t-1) * math.cos(5*t)


def crie_vetor(n_colunas):
    """ Função "crie_vetor". Recebe um único parâmetro e retorna um vetor nulo.

    Arguments:
        n_colunas {int} -- número de colunas ou dimensão do vetor

    Returns:
        list -- vetor nulo
    """    

    vetor = []
    for j in range(n_colunas):
        x = 0
        vetor.append(x)

    return vetor


def contorno(vetor, teste):
    """ Função "contorno". Recebe dois parâmetros e retorna um vetor com as condições de contorno aplicadas.

    Arguments:
        vetor {list} -- vetor nulo que será preenchido
        teste {string} -- variável que contém o tipo de teste que o usuário deseja executar

    Returns:
        list -- vetor com as condições de contorno do problema
    """    

    delta_x = 1 / n

    for i in range(len(vetor)):
        vetor[i] = u0(i * delta_x, teste)

    return vetor


main()
