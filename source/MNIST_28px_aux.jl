# ---------------------------------------------------------------------------
# Bibliotecas
# ---------------------------------------------------------------------------

# Instala bibliotecas externas. Por favor, descomente se necessario
#using Pkg
#pkg"add MLDatasets"
#pkg"add DataDeps"

using MLDatasets
using DataDeps

# ---------------------------------------------------------------------------
# Funcao para importar base de dados MNIST
# ---------------------------------------------------------------------------

#TODO - Comentar
function bring_me_the_MNIST(samples = 60000)
  #

  # Conjunto de treinamento train_x(28, 28, 60000) e train_y: (60000,)
  train_x, train_y = MNIST.traindata()

  # Conjunto de teste test_x: (28, 28, 10000) e test_y: (10000,)
  #test_x,  test_y  = MNIST.testdata()

  flatted_train_x = reshape(train_x, 784, 60000)

  flatted_train_x = Array(flatted_train_x[:, 1:samples]')

  return flatted_train_x, train_y[1:samples]
end

#TODO - Comentar
function show_me_the_MNIST(flatted_train_x, lin = 8, col = 4)
  #

  pixels = isqrt(size(flatted_train_x)[2])
  samples = size(flatted_train_x)[1]
  normal_train_x = reshape(flatted_train_x', (pixels, pixels, samples))
  images = MNIST.convert2image(normal_train_x[:, :, 1:(lin * col)])
  return mosaicview(images; fillvalue = 0, ncol = col, rowmajor = true)
end

#TODO - Comentar
function convert_to_image(flatted_train_x)
  #

  pixels = isqrt(size(flatted_train_x)[2])
  samples = size(flatted_train_x)[1]
  normal_train_x = reshape(flatted_train_x', (pixels, pixels, samples))
  images = MNIST.convert2image(normal_train_x)
  return images
end

# ------------------------------------------------------------------------------
# Funcao para localizar exemplos com alta e baixa variancias por componente
# ------------------------------------------------------------------------------
#TODO - Comentar
function fig_examples_relation_feature_space(A, B, U, V, n)
  #

  # Inicializa vetor de exemplos por componente
  examples=[]

  # Exemplos com maior e menor expressividade, componente a componente
  for i in 1:n

    # Indices dos maiores e menores fenomenos por componente e dataset
    min_B = argmin(V[:, i])
    max_B = argmax(V[:, i])
    min_A = argmin(U[:, i])
    max_A = argmax(U[:, i])

    # Calcula vetor de exemplos por componente
    push!(examples, [Array(B[max_B, :]), Array(B[min_B, :]), 
                     Array(A[min_A, :]), Array(A[max_A, :])])

  end

  return examples
  #return getindex(examples)
end

#TODO - Comentar
function examples_to_flatted_train(examples)
  #

  images = examples[1][1][:]
  for i in 1:length(examples)
    for j in 1:length(examples[i])
      images = hcat(images, examples[i][j][:])
    end
  end
  return images[:, 2:size(images)[2]]'
end

# ------------------------------------------------------------------------------
# Funcao para filtrar as amostras do MNIST
# ------------------------------------------------------------------------------

function filter_MNIST(data, label, filter)
  # Dado o dataset MNIST (amostras X pixels), sua classificacao e um valor de 
  # filtro, retorna um novo dataset apenas com as classificacoes, filtradas por
  # linhas.

  # Vetor de dados filtrados
  filtered_data = Array{N0f8}(undef, 0, size(data)[2])

  # Concatena as linhas com identificadores iguais ao filtro
  for i in 1:length(label)
    if label[i] == filter
      filtered_data = vcat(filtered_data, data[i,:]')
    end
  end

  # Retorna dados filtrados
  return Array(filtered_data)
end