# ---------------------------------------------------------------------------
# Bibliotecas
# ---------------------------------------------------------------------------

# Instala bibliotecas externas. Por favor, descomente se necessario
#using Pkg
#Pkg.add("CSV")
#Pkg.add("DataFrames")
#Pkg.add("ImageView")
using CSV
using DataFrames
using ImageView

# ---------------------------------------------------------------------------
# Funcao para importar base de dados
# ---------------------------------------------------------------------------
#TODO - Testar
#TODO - Comentar
function bring_me_the_MNIST_16px(filter = 10, database_path)
  #
  # UCI 1592x256             Classifications on [:,257:266]

  #
  df = CSV.read(string(database_path, "semeion.data"))

  #
  if       filter==0
    df =   df[df."1"   .> 0,:]
    elseif filter==1
      df = df[df."0"   .> 0,:]
    elseif filter==2
      df = df[df."0_1" .> 0,:]
    elseif filter==3
      df = df[df."0_2" .> 0,:]
    elseif filter==4
      df = df[df."0_3" .> 0,:]
    elseif filter==5
      df = df[df."0_4" .> 0,:]
    elseif filter==6
      df = df[df."0_5" .> 0,:]
    elseif filter==7
      df = df[df."0_6" .> 0,:]
    elseif filter==8
      df = df[df."0_7" .> 0,:]
    elseif filter==9
      df = df[df."0_8" .> 0,:]
  else
    df = df[:,:]
  end

  # 
  data = Matrix(df[:,1:256])

  # 
  label = Matrix(df[:,257:267])

  # 
  pixels = isqrt(size(data)[2])
  samples = size(data)[1]
  data_3D = reshape(data', (pixels, pixels, samples))
  images = colorview(Gray, data_3D)

  return data, label, images
end

#TODO - Avaliar se há algo para aproveitar e limpar
# ---------------------------------------------------------------------------
# Alguns testes temporarios
# ---------------------------------------------------------------------------

# df = CSV.read(string(database_path, "semeion.data"))
# data = Matrix(df[:,1:256])

# show(df[67,257:267])


# example_1 = data[4,:]
# pixels = isqrt(size(example_1)[1])
# image_1 = reshape(example_1,(pixels,pixels))
# imshow(image_1)

# example_2 = data[12,:]
# pixels = isqrt(size(example_2)[1])
# image_2 = reshape(example_2,(pixels,pixels))
# imshow(image_2)

# mosaicview(colorview(Gray,image_1),colorview(Gray,image_2); fillvalue=0, ncol=4, rowmajor=true)

# #images = colorview(Gray,data)