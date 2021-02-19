# ---------------------------------------------------------------------------
# Bibliotecas
# ---------------------------------------------------------------------------

include("GSVD_measures.jl")
include("GSVD_aux.jl")
include("GSVD_classifications.jl")
include("GSVD_plots.jl")
include("MORAES_studentXterms.jl")
include("MOVIELENS.jl")


# Instala bibliotecas externas. Por favor, descomente se necessario
# using Pkg
# pkg"add NMF"

using NMF # necessario para decomposicao NMF

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
# MODELO 1:  termos_alunos_questões          termos_conceitos_questões
#            ---[T+)----[Q)---        =        ---[Ct)----[Cq)---
#
# ou SVD(Q*T+)
#
#
# MODELO 2:  termos_alunos_questões          termos_conceitos_questões
#            ---[Tt)----[Q)---        =        ---[Ct)----[Cq)---
#
# ou SVD(Q*Tt)
#
#
# MODELO 3:  alunos_termos_conceitos          alunos_questões_conceitos
#            ---[T)----[Ct)---        =        ---[Q)----[Cq)---
#
# ou GSVD(T, Q)
#
#
# MODELO 4:  termos_alunos_questões          termos_conceitos_questões
#            ---[T+)----[Q)---        =        ---[Ct)----[Cq)---
#
# ou NMF(Q*T+)
#
#
# MODELO 5:  termos_alunos_questões          termos_conceitos_questões
#            ---[Tt)----[Q)---        =        ---[Ct)----[Cq)---
#
# ou NMF(Q*Tt)

T = [1  2  3
     4  5  6
     7  8  9
    10 11 12
    34 22 69
    63 46 78
    44 55 37
    34 37 87]
Q = [12 24 35
     43 57 68
     75 84 93
     59 88 17
     45 54 86
     66 88 99]

T, Q = bring_me_MOVIELENS(database_path)

function best_model(T, Q)
    # Calcula a diferença entre os modelos propostos

    # MODELO 1
    U, S, V = svd(Q * pinv(T))
    Ct = 
    Cq = 
    println("Modelo 1: ", norm((Q * pinv(T)) - (Cq * Ct))) 

    # MODELO 2
    U, S, V = svd(Q * T')
    Ct = 
    Cq = 
    println("Modelo 2: ", norm((Q * T')      - (Cq * Ct)))
    
    # MODELO 3
    U1, U2, E1, E2, X = normalized_gsvd(T, Q)
    Ct = U1
    Cq = U2
    println("Modelo 3: ", norm((Ct * T)      - (Cq * Q)))

    # MODELO 4
    H, W = nnmf(Q * pinv(T), 3)
    Ct = 
    Cq = 
    println("Modelo 4: ", norm((Q * pinv(T)) - (Cq * Ct))) 

    # MODELO 5
    H, W = nnmf(Q * T', 3)
    Ct = 
    Cq = 
    println("Modelo 5: ", norm((Q * T')      - (Cq * Ct)))
end

best_model(T, Q)