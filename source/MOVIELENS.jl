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

function bring_me_MOVIELENS(database_path)
    # Retorna os Dataframes e respectivas matrizes de filmes por genero e
    # filmes por usuario(1).
    #
    # (1) F. HARPER, J. KONSTAN. The MovieLens Datasets: History and Context.
    # ACM Transactions on Interactive Intelligent Systems (TiiS) 5, 4: 
    # 19:1–19:19.(2015) https://doi.org/10.1145/2827872
  
    # Importação de arquivos CSV em matrizes
    df_ratings = DataFrame(load(string(database_path, "ratings.csv"); 
                header_exists = true))
    df_movies = DataFrame(load(string(database_path, "movies.csv"); 
                header_exists = true))

    # Convertendo dataframes em formato wide
    
    # Filmes X Usuarios
    df_ratings_unstacked = unstack(df_ratings, :movieId, :userId, :rating);

    # Substituindo valores ausentes no Dataframe
    column_names = names(df_ratings_unstacked);
    for i in 1:length(column_names)
        df_ratings_unstacked[!, column_names[i]] = replace(
            df_ratings_unstacked[!, column_names[i]], missing => 0)
    end

    # Filmes X Genero

    # Inicializando matriz vazia da qual deve-se anexar os generos exclusivos
    # depois de dividi-los pelo separador '|'
    unique_genres = []

    # Selecionando generos unicos
    for i in df_movies.genres
        b = split(i, '|')
        for n in b
            if n ∉ unique_genres
                push!(unique_genres,n)
            end
        end
    end

    # Criando as colunas de antemão e as preenchendo com 0, evitando erro 
    # "ArgumentError: Cannot assign to non-existent column:" ao executar 
    # próximo código
    for x in unique_genres
        df_movies[!, Symbol(x)] .= 0
    end

    # Dividindo novamente os elementos da coluna de gêneros
    chris = split.(df_movies.genres, '|');

    # Substituindo o Dataframe com valores desejados (0 e 1)
    for x in unique_genres
        for ij in 1:length(df_movies.genres)
            if x ∉ chris[ij]
                df_movies[ij, Symbol(x)] = 0
            else
                df_movies[ij, Symbol(x)] = 1
            end
        end
    end

    # Retirando colunas de titulo e genero
    df_movies_unstacked = select(df_movies,           Not(:title))
    df_movies_unstacked = select(df_movies_unstacked, Not(:genres))

    # Selecionando apenas os filmes comuns entre ambos Dataframes, garantindo
    # mesmas dimensoes compartilhadas em ambas matrizes
    df_movies_unstacked = df_movies_unstacked[in.(df_movies_unstacked.movieId, Ref(df_ratings_unstacked[!, :movieId])), :]

    # Convertendo em matrizes
    ratings_matrix = df_ratings_unstacked[:,  2:size(df_ratings_unstacked)[2]]
    rating_matrix = convert(Matrix, ratings_matrix)
    movies_matrix = df_movies_unstacked[:,  2:size(df_movies_unstacked)[2]]
    movie_matrix = convert(Matrix, movies_matrix)

    # Transpondo para que a dimensão compartilhada esteja em colunas 
    rating_matrix = rating_matrix'
    movie_matrix = movie_matrix'

    # Retornando Dataframes e respectivas matrizes
    return df_ratings_unstacked, df_movies_unstacked, rating_matrix, movie_matrix
end