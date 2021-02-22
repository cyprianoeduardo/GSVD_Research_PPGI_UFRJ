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
function best_alg(T, Q, output_path, k1 = 0, k2 = 0)
    # Dados duas matrizes com uma dimensão em comum, verifica a precisão de 
    # Algoritmos com os Modelos propostos. 

    # Necessário para resolver o erro "ERROR: MethodError: no method 
    # matching eigen(::Adjoint{Float64,Array{Float64,2}})" em funcoes da 
    # biblioteca LinearAlgebra.
    # Ref.: https://github.com/JuliaLang/julia/issues/28714
    T = copy(T)
    Q = copy(Q)

    # Algoritmo 1
    println("Running Algorithm 1")
    U, S, V = svd(Q * pinv(T))
    if k1 == 0
        svd_knee(S, output_path)
        k1 = input("Please analyze the SVD knee graph and provide a value for K:")
    end
    Ct = sqrt(diagm(S[1:k1])) * V[:, 1:k1]'
    Cq = U[:, 1:k1] * sqrt(diagm(S[1:k1]))
    results = test_models(T, Q, Ct, Cq)

    # Algoritmo 2
    println("Running Algorithm 2")
    U, S, V = svd(Q * T')
    if k2 == 0
        svd_knee(S, output_path)
        k2 = input("Please analyze the SVD knee graph and provide a value for K:")
    end
    Ct = sqrt(diagm(S[1:k2])) * V[:, 1:k2]'
    Cq = U[:, 1:k2] * sqrt(diagm(S[1:k2]))
    results = vcat(results, test_models(T, Q, Ct, Cq))
    
    # REVIEW - Verificar como reduzir a dimensão
    # FIXME - Verificar dimensões de Ct e Cq, pois ocorre erro no primeiro modelo
    # Algoritmo 3
    println("Not Running Algorithm 3")
    # U, V, Q2, E1, E2, R = svd(T, Q)
    # X = R * Q2'
    # TQ = [U * E1; V * E2] * X
    # Ct = TQ[1:size(T)[1], :]
    # Cq = TQ[1:size(Q)[1], :]
    # println("Size T", size(T))
    # println("Size Q", size(Q))
    # println("Size Ct", size(Ct))
    # println("Size Cq", size(Cq))
    # results = vcat(results, test_models(T, Q, Ct, Cq))
    results = vcat(results, [Inf Inf Inf Inf])

    # Algoritmo 4
    println("Running Algorithm 4")
    try
        r = nnmf((Q * pinv(T)), k1; alg=:multmse, maxiter=30, tol=1.0e-4)
        W = r.W
        H = r.H
        Ct = H
        Cq = W
        results = vcat(results, test_models(T, Q, Ct, Cq))
    catch e
        println("Algorithm 4 didn't run because data isn't positive.")
        results = vcat(results, [Inf Inf Inf Inf])
    end

    # Algoritmo 5
    println("Running Algorithm 5")
    try
        r = nnmf((Q * T'), k2; alg=:multmse, maxiter=30, tol=1.0e-4)
        W = r.W
        H = r.H
        Ct = H
        Cq = W
        results = vcat(results, test_models(T, Q, Ct, Cq))
    catch e
        println("Algorithm 5 didn't run because data isn't positive.")
        results = vcat(results, [Inf Inf Inf Inf])
    end

    # Apresenta os resultados
    println(convert(DataFrame, results))
    println("The best is Algorithm ", argmin(results)[1], " in Model ", argmin(results)[2])

    return results
end

# Extraindo matrizes do Movielens
dfT, dfQ, T, Q = bring_me_MOVIELENS(database_path)

# Apresentando dimensoes
println("Size Users  X Movies:", size(T))
println("Size Genres X Movies:", size(Q))

# Definindo joelhos da curva
k1 = 0
k2 = 0

# Verificando resultados
println("- Raw data:")
results = best_alg(T, Q, output_path, k1, k2)

# ---------------------------------------------------------------------------
# Testando centralizacoes e normalizacoes das entradas
# ---------------------------------------------------------------------------

# Centralizando dados
T, T_means = centralizer(T)
Q, Q_means = centralizer(Q)

# Verificando resultados
println("- Centralized data:")
results = best_alg(T, Q, output_path, k1, k2)


# Centralizando stack
TQ = vcat(T, Q)
TQ, TQ_means = centralizer(TQ)
T = TQ[1:size(T)[1], :]
Q = TQ[1:size(Q)[1], :]

# Verificando resultados
println("- Centralized stack data:")
results = best_alg(T, Q, output_path, k1, k2)


# Normalizando dados
T, T_means, T_std = zscoretransform(T)
Q, Q_means, Q_std = zscoretransform(Q)

# Verificando resultados
println("- Standardized data:")
results = best_alg(T, Q, output_path, k1, k2)


# Normalizando stack
TQ = vcat(T, Q)
TQ, TQ_means, TQ_std = zscoretransform(TQ)
T = TQ[1:size(T)[1], :]
Q = TQ[1:size(Q)[1], :]

# Verificando resultados
println("- Standardized stack data:")
results = best_alg(T, Q, output_path, k1, k2)


# TODO - Testar com matrizes sintéticas