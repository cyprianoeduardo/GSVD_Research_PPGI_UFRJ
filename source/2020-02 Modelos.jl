# ---------------------------------------------------------------------------
# Bibliotecas
# ---------------------------------------------------------------------------

include("GSVD_measures.jl")
include("GSVD_aux.jl")
include("GSVD_classifications.jl")
include("GSVD_plots.jl")
include("SVD_plots.jl")
include("MORAES_studentXterms.jl")
include("MOVIELENS.jl")


# Instala bibliotecas externas. Por favor, descomente se necessario
# using Pkg
# # pkg"add NMF"
# pkg"add RCall"
# Pkg.build("RCall")

using NMF # necessario para decomposicao NMF
# using RCall # necessario para codigos na linguagem R

# ---------------------------------------------------------------------------
# Definindo caminhos para importação e exportação de arquivos
# ---------------------------------------------------------------------------

# Define local para exportação dos plots como imagens
output_path = string(pwd(), "/source/outputs/")

# Define local onde estão armazenados os dados
database_path = string(pwd(), "/source/databases/")

# ---------------------------------------------------------------------------
# Testes - Por favor, descomente conforme necessidade
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Padrão para rodar códigos na linguagem R
# ---------------------------------------------------------------------------

# @rput # Julia -> R

# R"""
# # code
# # goes
# # here
# """

# @rget # R -> Julia

# ---------------------------------------------------------------------------
# Reunião 12-02-2020
# ---------------------------------------------------------------------------

# INPUT:   alunos_termos        alunos_questões
#            ---[T)---            ---[Q)---
#
# OUTPUT:  termos_conceitos   conceitos_questões
#            ---[Ct)---           ---[Cq)---
#
# tal que |conceitos| << |alunos|
#
#
#
# Qual melhor modelo para alcançar o output?
#
# Modelo 1:  termos_alunos_questões          termos_conceitos_questões
#            ---[T+)----[Q)---        =        ---[Ct)----[Cq)---
#
# ou SVD(Q*T+) Algoritmo 1
# ou NMF(Q*T+) Algoritmo 4
#
#
# Modelo 2:  termos_alunos_questões          termos_conceitos_questões
#            ---[T')----[Q)---        =        ---[Ct)----[Cq)---
#
# ou SVD(Q*T') Algoritmo 2
# ou NMF(Q*T') Algoritmo 5
#
# Modelo 3:  alunos_termos_conceitos          alunos_questões_conceitos
#            ---[T)----[Ct)---        =        ---[Q)----[Cq+)---
#
# ou GSVD(T, Q) Algoritmo 3
#
# Modelo 4:  alunos_termos_conceitos          alunos_questões_conceitos
#            ---[T)----[Ct)---        =        ---[Q)----[Cq')---
#
# ou GSVD(T, Q) Algoritmo 3

# ---------------------------------------------------------------------------
# Funcao para teste de algoritmos com modelos propostos
# ---------------------------------------------------------------------------
function test_models(T, Q, Ct, Cq)
    # Retorna vetor com resultados dos modelos
    # println("Size T:", size(T))
    # println("Size Q:", size(Q))
    # println("Size pinv(T):", size(pinv(T)))
    # println("Size T':", size(T'))
    # println("Size Ct:", size(Ct))
    # println("Size Cq:", size(Cq))
    results = []
    push!(results, norm((Q  * pinv(T)) - (Cq       * Ct))) # Modelo 1
    push!(results, norm((Q  * T')      - (Cq       * Ct))) # Modelo 2
    push!(results, norm((Ct * T)       - (pinv(Cq) * Q ))) # Modelo 3
    push!(results, norm((Ct * T)       - (Cq'      * Q ))) # Modelo 4

    return results'
end

# ---------------------------------------------------------------------------
# Funcao para entrada de dados do usuario
# ---------------------------------------------------------------------------
function input(prompt::AbstractString="")::Int64
    print(prompt)
    return parse(Int64, chomp(readline()))
end

# ---------------------------------------------------------------------------
# Funcao para criar tabela de resultados de testes e mostrar melhor algoritmo
# ---------------------------------------------------------------------------
function best_alg(T, Q, output_path)
    # Calcula a diferença entre os modelos propostos

    # Algoritmo 1
    U, S, V = svd(Q * pinv(T))
    # println("Size U:", size(U))
    # println("Size diagm(S):", size(diagm(S)))
    # println("Size V:", size(V))
    # println("Size sqrt(S):", size(sqrt.(diagm(S))))
    svd_knee(S, output_path)
    k1 = input("Please analyze the SVD knee graph and provide a value for K:")
    Ct = sqrt.(diagm(S[1:k1])) * V[:, 1:k1]'
    Cq = U[:, 1:k1] * sqrt.(diagm(S[1:k1]))
    results = test_models(T, Q, Ct, Cq)

    # Algoritmo 2
    U, S, V = svd(Q * T')
    svd_knee(S, output_path)
    k2 = input("Please analyze the SVD knee graph and provide a value for K:")
    Ct = sqrt.(diagm(S[1:k2])) * V[:, 1:k2]'
    Cq = U[:, 1:k2] * diagm(S[1:k2])
    results = vcat(results, test_models(T, Q, Ct, Cq))
    
    # TODO - Verificar como escrever o algoritmo do GSVD
    # Algoritmo 3 # Não é prioridade
    # U1, U2, E1, E2, X = normalized_gsvd(T, Q)
    # Ct = U1 * E1 * X
    # Cq = U2 * E2
    # results = vcat(results, test_models(T, Q, Ct, Cq))

    # FIXME - pinv(T) gera números negativos! NNMF exige apenas valores positivos!
    # Algoritmo 4
    # H, W = nnmf(Q * pinv(T), k1)
    # Cq = H
    # Ct = W
    # results = vcat(results, test_models(T, Q, Ct, Cq))

    # FIXME - Erro eps() ?!
    # # Algoritmo 5
    # H, W = nnmf(Q * T', k2)
    # Cq = H
    # Ct = W
    # results = vcat(results, test_models(T, Q, Ct, Cq))

    # Apresenta os resultados
    println("O melhor é o Algoritmo ", argmin(results)[1], " no Modelo ", argmin(results)[2])
    println(convert(DataFrame, results))

    return results
end

# Extraindo matrizes do Movielens
dfT, dfQ, T, Q = bring_me_MOVIELENS(database_path)

# Apresentando dimensões
println("Usuários X Filmes:", size(T))
println("Gêneros X Filmes:", size(Q))

# Verificando resultados
results = best_alg(T, Q, output_path)

# TODO - Verificar onde alocar a matriz S. Calcular a raiz e inserir em Ct e Cq?
# TODO - Testar centralizações e normalizações das entradas
# TODO - Testar com matrizes sintéticas