# ---------------------------------------------------------------------------
# Bibliotecas
# ---------------------------------------------------------------------------

include("GSVD_measures.jl")
include("GSVD_aux.jl")
include("GSVD_classifications.jl")
include("GSVD_plots.jl")
include("MORAES_studentXterms.jl")
include("MNIST_28px_aux.jl")

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
# Testes com classificador
# ---------------------------------------------------------------------------

#test_classifications(5, 5, 4) #AltaxAlta # OK!
#test_classifications(3, 3, 4) #GordaxGorda # FIXME - Erro na classificação de matrizes gordas 
#test_classifications(5, 3, 4) #AltaxGorda # FIXME - Erro na classificação de matrizes gordas 
#test_classifications(3, 5, 4) #GordaxAlta # FIXME - Erro na classificação de matrizes gordas 

# ---------------------------------------------------------------------------
# Testes com matrizes simples
# ---------------------------------------------------------------------------

# A = [1  2  3
#      4  5  6
#      7  8  9
#     10 11 12
#     34 22 69
#     63 46 78
#     44 55 37
#     34 37 87]
# B = [12 24 35
#      43 57 68
#      75 84 93
#      59 88 17
#      45 54 86
#      66 88 99]

# ---------------------------------------------------------------------------
# Testes com matrizes randomicas
# ---------------------------------------------------------------------------

# A = rand(6000, 300)
# B = rand(6000, 300)

# ---------------------------------------------------------------------------
# Testes com Datasets Laura
# ---------------------------------------------------------------------------

# A, B = bring_me_studentsXterms(database_path)

# ---------------------------------------------------------------------------
# Testes com Dataset MNIST
# ---------------------------------------------------------------------------

data, label = bring_me_the_MNIST()

# Filtrando MNIST por numeros
A = filter_data(data, label, 3)
B = filter_data(data, label, 8)

# ---------------------------------------------------------------------------
# Testes com GSVD
# ---------------------------------------------------------------------------

# GSVD Normalizado
U, V, D1, D2, X = normalized_gsvd(A, B)

# # Imprime dimensoes
# println("Size A: ", size(A))
# println("Size B: ", size(B))
# println("Size U: ", size(U))
# println("Size V: ", size(V))
# println("Size D1: ", size(D1))
# println("Size D2: ", size(D2))
# println("Size X: ", size(X))

# ---------------------------------------------------------------------------
# Testes com metricas do GSVD
# ---------------------------------------------------------------------------

# Plota metricas
show_relation_measures(D1, D2, output_path)

# # Valores singulares e ordem dos componentes pos ordenamento
# #alphas, betas, gsv_sorted_idx = alphas_betas(D1, D2, true)

# # Calcula metricas
# theta = antisymmetric_angular_distance(alphas./betas)
# P1, P2 = generalized_fractions_eigenexpression(alphas, betas)
# gsv = alphas./betas

# # Imprime dimensoes
# println("Size alphas: ", size(alphas))
# println("Size betas: ", size(betas))
# println("Size gsv_sorted_idx: ", size(gsv_sorted_idx))
# println("Size theta: ", size(theta))
# println("Size P1: ", size(P1))
# println("Size P2: ", size(P2))
# println("Size gsv: ", size(gsv))

# ---------------------------------------------------------------------------
# Testes com plots do MNIST
# ---------------------------------------------------------------------------

# shared_dim_size = size(X)[2]

# X = fig_examples_relation_feature_space(A, B, U, V, shared_dim_size)

# Y = examples_to_flatted_train(X)

# save(string(output_path, "examples.png"), show_me_the_MNIST(Y', length(X), 4))