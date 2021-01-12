# ---------------------------------------------------------------------------
# Biblioteca Intel MKL BLAS
# ---------------------------------------------------------------------------

# Devido a problemas anteriores com os resultados do GSVD, decorrentes de 
# versoes antigas do LAPACK no OpenBLAS, estamos utilizando o Intel MKL BLAS
# (<https://github.com/JuliaComputing/MKL.jl>)
# Atentar que a versão 1.5 apresenta erro ao compilar o MKL. Portanto, 
# estamos utilizando o Julia versão 1.4.

# Instala o Intel MKL BLAS. Por favor, descomente se necessario.
#using Pkg
#pkg"add MKL"
#pkg"build MKL"

# Testa instalação do MKL
using LinearAlgebra
println(BLAS.vendor())

# ---------------------------------------------------------------------------
# Bibliotecas
# ---------------------------------------------------------------------------

using LinearAlgebra #necessario para funcoes de Algebra Linear
using Random # necessario para funcoes sort

# ---------------------------------------------------------------------------
# Funcao para dividir conjunto de dados em treino e teste
# ---------------------------------------------------------------------------

function partitionTrainTest(data, at = 0.7)
    # Dado um conjunto de dados e um percentual, divide um conjunto de dados
    # entre teste e treino
    n = size(data)[1]
    idx = shuffle(1:n)
    train_idx = view(idx, 1:floor(Int, at * n))
    test_idx = view(idx, (floor(Int, at * n) + 1):n)
    return data[train_idx, :], data[test_idx, :]
end

# ---------------------------------------------------------------------------
# Funcao para normalizar matrizes da fatoração GSVD
# ---------------------------------------------------------------------------
#TODO - Finalizar comentários
function normalized_gsvd(A, B)
    # Dados duas matrizes A e B, executa o GSVD(A, B), normalizando os 
    # Componentes para que formem uma base.
    #
    # Adaptado de "GSVDAnalysis.m"(1).
    #
    # 1. FERGUSON, A. Exploring the Yeast Genome with Generalized Singular 
    # Value Decomposition. 2008. Princeton University, 2008. 
    # DOI: 10.1145/2480362.2480536.

    # Calcula o GSVD de A e B, onde :
    # X = (R * Q')'
    # A = U1 * E1 * X'
    # B = U2 * E2 * X'
    U1, U2, Q, E1, E2, R = svd(A, B)

    X = (R * Q')'
    
    X = X'     # " GSVD = U * E * X' but we want to work with U * E * X "
    
    # " X is now genelets x arrays [rows x cols] " 
    # " E(1,2) is now arraylets x genelets "
    # " U(1,2) is now genes x arraylets "
    
    X_rows = size(X)[1]
    X_cols = size(X)[2]
    
    # " Next; normalize the genelets so that we have a set of equal-length() 
    # basis vectors. This requires scaling the values of E(1,2) so that 
    # we maintain M(1,2) = U(1,2) * E(1,2) * X "

    normedX  = zeros(X_rows, X_cols) 
    normedE1 = zeros(size(E1)[1], size(E1)[2])
    normedE2 = zeros(size(E2)[1], size(E2)[2])

    # Localiza onde inicia a diagonal em D2
    k = argmax(E2[1, :]) - 1
    # Extrai os valores singulares das respectivas diagonais
    betas = Array(diag(E2, k))
    # Calcula vetores de valores nulos. As quantidades se baseiam na diferenca
    # entre o tamanho da dimensao compartilhada de A e B e o quantidade de 
    # valores singulares presentes nas diagonais de D1 e D2.
    zeros_beta = zeros(size(E2)[2] - length(betas))
    # Atualiza valores singulares com os valores nulos
    betas  = vcat(zeros_beta, betas)
    # Reconstroi E2
    E2 = diagm(betas)
    
    
    for i in 1:size(X)[1]
        len = norm(X[i, 1:X_cols])
    
        normedX[i, 1:X_cols] = X[i,1:X_cols] ./ len
        
        normedE1[i, i] = E1[i, i] * len
        normedE2[i, i] = E2[i, i] * len
    end
    
    X  = normedX
    E1 = normedE1
    E2 = normedE2

    return U1, U2, E1, E2, X
end
