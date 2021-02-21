# ---------------------------------------------------------------------------
# Bibliotecas
# ---------------------------------------------------------------------------

# Instala bibliotecas externas. Por favor, descomente se necessario
#using Pkg
#pkg"add Plots"
#pkg"add Images"
#pkg"add ORCA"
#pkg"add PlotlyBase"

using Plots # necessario para plot de graficos
using Images # necessario para manipulação de imagens
using ORCA # necessario para salvar plots em .png 


# ---------------------------------------------------------------------------
# Funcao para visualizar o "joelho" da curva
# ---------------------------------------------------------------------------

function svd_knee(S, output_path, log = false)
    # Dada o vetor de valores singulares do SVD, plota uma representacao
    # grafica da relevância dos principais componentes. É possível 
    # representa-lo numa escala logaritmica.

    # Inicializa motor grafico de plotagem
    plotly()

    # Define se usara escala logaritmica
    if log == true
        S = log.(S)
    end

    # Plot da relevancia
    p1 = plot(range(1, length=length(S)), S,
      title = "Joelho do SVD",
      xlabel = "K componentes principais",
      ylabel = "Relevância",
      xrotation = 60,
      leg = false,
      titlefontsize = 7,
      title_location = :center,
      tickfontsize = 6,
      legendfontsize = 6,
      guidefontsize = 6,
      legendtitlefontsize = 6
    )

    plot(p1, dpi = 600)
    png(string(output_path, "svd_knee.png"))
end