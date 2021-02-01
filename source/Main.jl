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

# Extraindo o MNIST
data, label = bring_me_the_MNIST()

number1 = 3
number2 = 8

# Filtrando MNIST por numeros
A = filter_MNIST(data, label, number1)
B = filter_MNIST(data, label, number2)

# Centralizando datasets distintamente
A, A_means = centralizer(A)
B, B_means = centralizer(B)

# ---------------------------------------------------------------------------
# Testes com GSVD
# ---------------------------------------------------------------------------

# GSVD Normalizado
U, V, D1, D2, X = normalized_gsvd(A, B)

# ---------------------------------------------------------------------------
# Testes com metricas do GSVD
# ---------------------------------------------------------------------------

# Plota metricas
show_relation_measures(D1, D2, output_path)

# Valores singulares e ordem dos componentes pos ordenamento
alphas, betas, gsv_sorted_idx = alphas_betas(D1, D2, true)

# Calcula metricas
theta = antisymmetric_angular_distance(alphas./betas)
P1, P2 = generalized_fractions_eigenexpression(alphas, betas)
gsv = alphas./betas

# ---------------------------------------------------------------------------
# Testes com plots do MNIST
# ---------------------------------------------------------------------------

# # Levantando o tamanho da dimensão compartilhada
# shared_dim_size = size(X)[2]

# # Descentralizando A e B, para correto plot das imagens de exemplos.
# A = filter_MNIST(data, label, number1)
# B = filter_MNIST(data, label, number2)

# # Extraindo os exemplos de amostras com maiores e menores variancias nos
# # componentes, tanto para o dataset A quanto para o B
# X2 = fig_examples_relation_feature_space(A, B, U, V, shared_dim_size)

# # Reajustando exemplos no formato amostra por pixels
# Y = examples_to_flatted_train(X2)

# # Plotando exemplos
# save(string(output_path, "A_B_examples.png"), show_me_the_MNIST(Y, shared_dim_size, 4))

# ---------------------------------------------------------------------------
# Escolhendo e plotando pinceis e respectivos exemplos
# ---------------------------------------------------------------------------

# # Escolhendo os melhores valores singulares 
# idx = best_sv_idx(1, alpha, betas, theta, P1, P2, pi / 32)

# Escolhendo componente com distancia angular proxima a zero.
idx_half = findfirst(x->(isapprox(x, 0; atol = pi / 256)), theta)

# Descentralizando A e B, para correto plot das imagens de exemplos.
A = filter_MNIST(data, label, number1)
B = filter_MNIST(data, label, number2)

# Plotando pincel com distancia angular proxima a zero e seu exemplo
save(string(output_path, "brush_A_B_angular_",    idx_half,      "_component.png"),          map(clamp01nan, convert_to_image(X)[:, :, idx_half]))
save(string(output_path, "brush_A_B_angular_",    idx_half,      "_component exampleA.png"),                 convert_to_image(A)[:, :, argmax(U[:, idx_half])])
save(string(output_path, "brush_A_B_angular_",    idx_half,      "_component exampleB.png"),                 convert_to_image(B)[:, :, argmax(V[:, idx_half])])

# Plotando pincel com distancia angular maxima para A e seu exemplo
save(string(output_path, "brush_A_max_angular_",  argmax(theta), "_component.png"),          map(clamp01nan, convert_to_image(X)[:, :, argmax(theta)]))
save(string(output_path, "brush_A_max_angular_",  argmax(theta), "_component example.png"),                  convert_to_image(A)[:, :, argmax(U[:, argmax(theta)])])

# Plotando pincel com distancia angular maxima para B e seu exemplo
save(string(output_path, "brush_B_max_angular_",  argmin(theta), "_component.png"),          map(clamp01nan, convert_to_image(X)[:, :, argmin(theta)]))
save(string(output_path, "brush_B_max_angular_",  argmin(theta), "_component example.png"),                  convert_to_image(B)[:, :, argmax(V[:, argmin(theta)])])

# Plotando pincel com fracao de autoexpressao maxima para A e seu exemplo
save(string(output_path, "brush_A_max_fraction_", argmax(P1),    "_component.png"),          map(clamp01nan, convert_to_image(X)[:, :, argmax(P1)]))
save(string(output_path, "brush_A_max_fraction_", argmax(P1),    "_component example.png"),                  convert_to_image(A)[:, :, argmax(U[:, argmax(P1)])])

# Plotando pincel com fracao de autoexpressao maxima para B e seu exemplo
save(string(output_path, "brush_B_max_fraction_", argmax(P2),    "_component.png"),          map(clamp01nan, convert_to_image(X)[:, :, argmax(P2)]))
save(string(output_path, "brush_B_max_fraction_", argmax(P2),    "_component example.png"),                  convert_to_image(B)[:, :, argmax(V[:, argmax(P2)])])