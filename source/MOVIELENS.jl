# ---------------------------------------------------------------------------
# Bibliotecas
# ---------------------------------------------------------------------------

# Instala bibliotecas externas. Por favor, descomente se necessario
# using Pkg
# pkg"add CSVFiles"
# pkg"add DataFrames"

using CSVFiles
using DataFrames

# ---------------------------------------------------------------------------
# Funcao para importar base de dados
# ---------------------------------------------------------------------------

database_path = string(pwd(), "/source/databases/")

# function bring_me_MOVIELENS(database_path)
#     # Retorna as matrizes de filmes por gênero e filmes por usuario(1).
#     #
#     #
#     # (1) F. HARPER, J. KONSTAN. The MovieLens Datasets: History and Context.
#     # ACM Transactions on Interactive Intelligent Systems (TiiS) 5, 4: 
#     # 19:1–19:19.(2015) https://doi.org/10.1145/2827872
  
    # Importação de arquivos CSV em matrizes
    df_ratings = DataFrame(load(string(database_path, "ratings.csv"); header_exists = true))
    df_movies = DataFrame(load(string(database_path, "movies.csv"); header_exists = true))

    # Convertendo dataframes em formato wide
    
    # Filmes X Usuarios
    df_ratings_unstacked = unstack(df_ratings, :movieId, :userId, :rating);
    # replace the missing values in the dataframe
    column_names_ru = names(df_ratings_unstacked);
    for i in 1:length(column_names_ru)
        df_ratings_unstacked[!, column_names_ru[i]] = replace(df_ratings_unstacked[!, column_names_ru[i]], missing => 0) # If we want to replace the missing values with any other number, then we use: missing => NaN
    end

    # Filmes X Genero
    # Criamos uma matriz vazia da qual devemos anexar em gêneros exclusivos depois de dividi-los pelo '|' separador
    l = []
    for i in df_movies.genres
        b = split(i, '|')
        for n in b
            if n ∉ l
                push!(l,n)
            end
        end
    end
    # Agora temos nossa variedade de gêneros únicos

    # Nós criamos as colunas de antemão e as preenchemos com 0's dos quais, se não o fizéssemos, terminaríamos com um erro "ArgumentError: Cannot assign to non-existent column:" ao executar nosso próximo código
    for x in l
        df_movies[!, Symbol(x)] .= 0
    end
    # Nós novamente dividimos os elementos da coluna de gêneros
    chris = split.(df_movies.genres, '|');
    # Agora, substituímos o dataframe com nossos valores desejados
    for x in l
        for ij in 1:length(df_movies.genres)
            if x ∉ chris[ij]
                df_movies[ij, Symbol(x)] = 0
            else
                df_movies[ij, Symbol(x)] = 1
            end
        end
    end
    
    df_movies_unstacked = select(df_movies, Not(:title), Not(:genres))
  
#     return Matrix{Float64}(df1), Matrix{Float64}(df2)
# end

#first(df_ratings_unstacked)
first(df_movies_unstacked)