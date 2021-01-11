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

function bring_me_studentsXterms(database_path)
    # Retorna as matrizes de estudantes e termos(1).
    #
    # 1. MORAES, L. PEDREIRA, C. Clustering Introductory Computer Science
    # Exercises Using Topic Modeling Techniques. 2020. Systems and Computing 
    # Engineering (COPPE-PESC), Universidade Federal do Rio de Janeiro (UFRJ)
    # , 2020.
  
    # Importação de arquivos CSV em matrizes
    df1 = DataFrame(load(string(database_path, "student_performance.csv"); header_exists = false))
    df2 = DataFrame(load(string(database_path, "terms.csv");               header_exists = false))
  
    return Matrix{Float64}(df1), Matrix{Float64}(df2)
end