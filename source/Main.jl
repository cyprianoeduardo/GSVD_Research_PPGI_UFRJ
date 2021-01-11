# ---------------------------------------------------------------------------
# Bibliotecas
# ---------------------------------------------------------------------------

include("GSVD_measures.jl")
include("GSVD_aux.jl")
include("GSVD_classifications.jl")
include("GSVD_plots.jl")
include("MNIST_28px_aux.jl")

# ---------------------------------------------------------------------------
# Definindo caminhos para importação e exportação de arquivos
# ---------------------------------------------------------------------------

# Define local para exportação dos plots como imagens
output_path = string(pwd(), "/source/outputs/")

# Define local onde estão armazenados os dados
database_path = string(pwd(), "/source/databases/")

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
# Testes com matrizes randômicas
# ---------------------------------------------------------------------------

# A = rand(6000, 300)
# B = rand(6000, 300)

# ---------------------------------------------------------------------------
# Testes com Datasets Laura
# ---------------------------------------------------------------------------

# Instala bibliotecas externas. Por favor, descomente se necessario
# using Pkg
# pkg"add CSVFiles"
# pkg"add DataFrames"

# using CSVFiles
# using DataFrames

# df1 = DataFrame(load("student_performance.csv"; header_exists = false))
# df2 = DataFrame(load("terms.csv"; header_exists = false))
# A = Matrix{Float64}(df1)
# B = Matrix{Float64}(df2)

# ---------------------------------------------------------------------------
# Testes com Dataset MNIST
# ---------------------------------------------------------------------------

A,a = bring_me_the_MNIST()

A_train = filter_data(A, a, 3)
B_test = filter_data(A, a, 8)

A = A_train
B = B_test





show_relation_measures(D1, D2, output_path)

# # Valores singulares e ordem dos componentes pos ordenamento
# #alphas, betas, gsv_sorted_idx1 = alphas_betas(D1, D2, true)
# alphas = diag(D1) #REVIEW - ATUALIZAR NO GSVD_measures
# betas  = diag(D2) #REVIEW - ATUALIZAR NO GSVD_measures

# # Calcula metricas
# theta = antisymmetric_angular_distance(alphas./betas)
# P1, P2 = generalized_fractions_eigenexpression(alphas, betas)
# gsv = alphas./betas

# Imprime dimensoes
println("Size A: ", size(A))
println("Size B: ", size(B))
println("Size U: ", size(U))
println("Size V: ", size(V))
println("Size Q: ", size(Q))
println("Size D1: ", size(D1))
println("Size D2: ", size(D2))
println("Size R: ", size(R))
println("Size X: ", size((R * Q')))
# println("Size alphas: ", size(alphas))
# println("Size betas: ", size(betas))

# Plota metricas
show_relation_measures(D1, D2)




X = (R * Q')' 
X = X'
shared_dim_size = size(X)[2]

X = fig_examples_relation_feature_space(A, B, U, V, shared_dim_size)

Y = examples_to_flatted_train(X)

save("examples.png",show_me_the_MNIST(Y',length(X),4))