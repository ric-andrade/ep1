<!-- Este arquivo foi escrito em markdown, para uma melhor visualização utilize algum programa ou site que interprete a extensão .md -->

# Primeiro exercício programa - Disciplina MAP3121 - Turmas 5 e 6

---

## Alunos:

* Guilherme Fernandes Alves - 10774361
* Ricardo Aguiar de Andrade - 10774674

---

## Linguagem de programação utilizada:

> **Python 3**

---

## Bibliotecas utilizadas:

1. math
2. numpy
3. matplotlib

Todas em sua última versão.

---

## Considerações iniciais

1. O comprimento da barra foi normalizado para 1.
2. O intervalo total de tempo foi fixado em T=1s (porém o programa foi arquitetado para receber qualquer valor positivo de T, sendo necessário alterar o código fonte para habilitar o input ao usuário)

---
## Instruções de uso:

Ao executar o programa, a primeira informação pedida ao usuário será o método a ser utilizado. Neste exercício programa haverão três métodos:

* Euler explícito: *d*
* Euler implícito: *e*
* Crank-Nicolson: *k*

> Portanto, o primeiro parâmetro que o usuário deve fornecer é do tipo *string*.

O segundo parâmetro é o número **inteiro** *N* que representa a quantidade de partições do comprimento da barra estudada.

> Portanto, o segundo parâmetro que o usuário deve fornecer é do tipo *int*.

No caso do usuário escolher o método de Euler explícito (*d*) também será pedido o parâmetro $\lambda$. Este parâmetro pode ser qualquer número real positivo diferente de $0$. 

> Portanto, o terceiro parâmetro que o usuário deve fornecer é do tipo *float*, **condicionado ao caso do método escolhido ser Euler explícito**.

Vale notar que este método se comporta de maneira estável somente se $\lambda \leq 0.5$. Caso isso não ocorra poderão ocorrer situações em que um erro *nan* aparecerá na tela, indicando complicações nos cálculos realizados. Se isso acontecer, serão exibidas respostas degeneradas. 

Em seguida será pedido ao usuário qual teste ele deseja executar. Temos três possíveis testes:

* Teste a: *a*
* Teste b: *b*
* Teste c: *c*

> Portanto, o quarto parâmetro que o usuário deve fornecer é do tipo *string*.

Cada teste define uma solução exata, função fonte e condições de contorno, não necessariamente iguais. Para saber mais detalhes de cada teste consulte o relatório deste exercício programa. 

Por fim, a última informação que o usuário deve fornecer é se ele quer ou não visualizar os gráficos:

* Se sim: *y*
* Se não: *n*

> Portanto, o quinto parâmetro que o usuário deve fornecer é do tipo *string*.

Caso o usuário opte por ver os gráficos, estes serão exibidos em intervalos de 0.1s até 1.0s, ou seja, serão apresentados na tela 10 gráficos de temperatura $\times$ posição. A cada gráfico exibido também será mostrado o gráfico de erro $\times$ posição para aquele instante de tempo, totalizando 20 gráficos para cada rodada no código.

Em resumo, o primeiro gráfico se referi-rá ao instante 0.1s; ao fechar este gráfico será mostrado o gráfico dos erros para este instante; para passar para o instante de tempo 0.2s basta fechar este gráfico exibido e assim sucessivamente.

Ao fim da exibição dos gráficos serão mostrados os vetores de temperaturas aproximadas, temperaturas exatas e erros, bem como o erro máximo para o último instante de tempo, isto é t=1s. 

### Fim