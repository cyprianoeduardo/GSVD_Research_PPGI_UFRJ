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
# Reunião 23-02-2020
# ---------------------------------------------------------------------------

# INPUT:
#
# A = Usuarios X Filmes
# B = Generos X Filmes
#
# OUTPUT:
#
# Ca = Conceitos X Usuarios
# Cb = Conceitos X Generos
#
# MODELOS:
#                    G X F   F X U      G X C      C X U
# Modelo 1: MIN norm((B    * pinv(A)) - (pinv(Cb)*  Ca))
#
#                    C X U   U X F      C X G      G X F
# Modelo 2: MIN norm(Ca    * A)       - (Cb      * B )))

# ---------------------------------------------------------------------------
# Funcao para teste de algoritmos com modelos propostos
# ---------------------------------------------------------------------------
function test_models(A, B, Ca, Cb)
    # Retorna vetor com resultados dos modelos

    results = []
    push!(results, norm((B  * pinv(A)) - (pinv(Cb) * Ca))) # Modelo 1
    push!(results, norm((Ca * A)       - (Cb       * B ))) # Modelo 2

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
function best_alg(A, B, output_path, k = 0, nnegative = false, qt_models = 2)
    # Dados duas matrizes com uma dimensão em comum, verifica a precisão de 
    # Algoritmos com os Modelos propostos. 

    # Necessário para resolver o erro "ERROR: MethodError: no method 
    # matching function(::Adjoint{Float64,Array{Float64,2}})" em funcoes da 
    # biblioteca LinearAlgebra.
    # Ref.: https://github.com/JuliaLang/julia/issues/28714
    A = copy(A)
    B = copy(B)

    # Algoritmo 1
    println("Running Algorithm 1")
    U, S, V = svd(B * pinv(A))
    if k == 0
        svd_knee(S, output_path)
        k = input("Please analyze the SVD knee graph and provide a value for K:")
    end
    # Ca = sqrt(diagm(S[1:k])) * V[:, 1:k]'
    # Cb = U[:, 1:k] * sqrt(diagm(S[1:k]))
    Ca = V[:, 1:k]'
    Cb = pinv(U[:, 1:k])
    results = test_models(A, B, Ca, Cb)

    # Algoritmo 2
    println("Running Algorithm 2")
    try
        if nnegative == true
            r = nnmf(abs.(B * pinv(A)), k; alg=:projals, maxiter=30, tol=1.0e-4)
        else
            r = nnmf(    (B * pinv(A)), k; alg=:projals, maxiter=30, tol=1.0e-4)
        end
        W = r.W
        H = r.H
        Ca = H
        Cb = pinv(W)
        results = vcat(results, test_models(A, B, Ca, Cb))
    catch e
        println("Algorithm 2 didn't run because data isn't positive.")
        results = vcat(results, fill(Inf, (1, qt_models)))
    end

    # Algoritmo 3
    # REVIEW - Verificar como reduzir a dimensão
    # FIXME - Verificar dimensões de Ca e Cb, pois ocorre erro no primeiro modelo
    # Algoritmo 3
    println("Not Running Algorithm 3")
    # U, V, Q, E1, E2, R = svd(A, B)
    # X = R * Q'
    # AB = [U * E1; V * E2] * X
    # Ca = AB[1:size(A)[1], :]
    # Cb = AB[1:size(B)[1], :]
    # println("Size A", size(A))
    # println("Size B", size(B))
    # println("Size Ca", size(Ca))
    # println("Size Cb", size(Cb))
    # results = vcat(results, test_models(A, B, Ca, Cb))
    results = vcat(results, fill(Inf, (1, qt_models)))

    # Algoritmo 4 - 1.1
    println("Running Algorithm 1.1 pinv(A) e U'")
    U, S, V = svd(B * pinv(A))
    if k == 0
        svd_knee(S, output_path)
        k = input("Please analyze the SVD knee graph and provide a value for K:")
    end
    Ca = V[:, 1:k]'
    Cb = U[:, 1:k]'
    results = vcat(results, test_models(A, B, Ca, Cb))

    # Algoritmo 5 - 1.2
    println("Running Algorithm 1.2 A' e pinv(U)")
    U, S, V = svd(B * A')
    if k == 0
        svd_knee(S, output_path)
        k = input("Please analyze the SVD knee graph and provide a value for K:")
    end
    Ca = V[:, 1:k]'
    Cb = pinv(U[:, 1:k])
    results = vcat(results, test_models(A, B, Ca, Cb))

    # Algoritmo 6 - 1.3
    println("Running Algorithm 1.3 A' e U'")
    U, S, V = svd(B * A')
    if k == 0
        svd_knee(S, output_path)
        k = input("Please analyze the SVD knee graph and provide a value for K:")
    end
    Ca = V[:, 1:k]'
    Cb = U[:, 1:k]'
    results = vcat(results, test_models(A, B, Ca, Cb))

    # Algoritmo 7 - 2.1
    println("Running Algorithm 2.1 pinv(A) e W'")
    try
        if nnegative == true
            r = nnmf(abs.(B * pinv(A)), k; alg=:projals, maxiter=30, tol=1.0e-4)
        else
            r = nnmf(    (B * pinv(A)), k; alg=:projals, maxiter=30, tol=1.0e-4)
        end
        W = r.W
        H = r.H
        Ca = H
        Cb = W'
        results = vcat(results, test_models(A, B, Ca, Cb))
    catch e
        println("Algorithm 2.1 didn't run because data isn't positive.")
        results = vcat(results, fill(Inf, (1, qt_models)))
    end

    # Algoritmo 8 - 2.2
    println("Running Algorithm 2.2 A' e pinv(W)")
    try
        if nnegative == true
            r = nnmf(abs.(B * A'), k; alg=:projals, maxiter=30, tol=1.0e-4)
        else
            r = nnmf(    (B * A'), k; alg=:projals, maxiter=30, tol=1.0e-4)
        end
        W = r.W
        H = r.H
        Ca = H
        Cb = pinv(W)
        results = vcat(results, test_models(A, B, Ca, Cb))
    catch e
        println("Algorithm 2.2 didn't run because data isn't positive.")
        results = vcat(results, fill(Inf, (1, qt_models)))
    end

    # Algoritmo 9 - 2.3
    println("Running Algorithm 2.3 A' e W'")
    try
        if nnegative == true
            r = nnmf(abs.(B * A'), k; alg=:projals, maxiter=30, tol=1.0e-4)
        else
            r = nnmf(    (B * A'), k; alg=:projals, maxiter=30, tol=1.0e-4)
        end
        W = r.W
        H = r.H
        Ca = H
        Cb = W'
        results = vcat(results, test_models(A, B, Ca, Cb))
    catch e
        println("Algorithm 2.3 didn't run because data isn't positive.")
        results = vcat(results, fill(Inf, (1, qt_models)))
    end

    # Apresenta os resultados
    println(convert(DataFrame, results))
    println("The best is Algorithm ", argmin(results)[1], " in Model ", argmin(results)[2])

    return results
end

# Extraindo matrizes do Movielens
dfA, dfB, A, B = bring_me_MOVIELENS(database_path)

# Apresentando dimensoes
println("Size Users  X Movies:", size(A))
println("Size Genres X Movies:", size(B))

# Definindo joelho da curva
k = 20
nnegative = true

# Verificando resultados
println("- Raw data:")
results = best_alg(A, B, output_path, k, nnegative)

# ---------------------------------------------------------------------------
# Testando centralizacoes e normalizacoes das entradas
# ---------------------------------------------------------------------------

# Centralizando dados
A, A_means = centralizer(A)
B, B_means = centralizer(B)

# Verificando resultados
println("- Centralized data:")
results = best_alg(A, B, output_path, k, nnegative)


# Centralizando stack
AB = vcat(A, B)
AB, AB_means = centralizer(AB)
A = AB[           1:size(A )[1], :]
B = AB[size(A)[1]+1:size(AB)[1], :]

# Verificando resultados
println("- Centralized stack data:")
results = best_alg(A, B, output_path, k, nnegative)


# Normalizando dados
A, A_means, A_std = zscoretransform(A)
B, B_means, B_std = zscoretransform(B)

# Verificando resultados
println("- Standardized data:")
results = best_alg(A, B, output_path, k, nnegative)


# Normalizando stack
AB = vcat(A, B)
AB, AB_means, AB_std = zscoretransform(AB)
A = AB[           1:size(A )[1], :]
B = AB[size(A)[1]+1:size(AB)[1], :]

# Verificando resultados
println("- Standardized stack data:")
results = best_alg(A, B, output_path, k, nnegative)


# TODO - Testar com matrizes sintéticas