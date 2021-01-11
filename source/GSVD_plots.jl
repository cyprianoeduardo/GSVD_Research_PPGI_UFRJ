# ---------------------------------------------------------------------------
# Bibliotecas
# ---------------------------------------------------------------------------

# Instala bibliotecas externas. Por favor, descomente se necessario
#using Pkg
#pkg"add Plots"
#pkg"add Images"
#pkg"add ImageCore"
#pkg"add ImageShow"
#pkg"add ORCA"
#pkg"add PlotlyBase"

using Plots # necessario para plot de graficos
using Images # necessario para manipulação de imagens
using ImageCore # necessario para a funcao convert2image
using ImageShow # necessario para imagens em mosaico de grandes densidades
using ORCA # necessario para salvar plots em .png 


# ---------------------------------------------------------------------------
# Funcao para visualizacoes para as funcoes Distancia Angular Antisimetrica e
# Fracoes Generalizadas de Autoexpressao
# ---------------------------------------------------------------------------

function show_relation_measures(D1, D2, output_path)
    # Dadas as matrizes D1 e D2 da fatoracao GSVD(A,B), plota as metricas da 
    # Distancia Angular Antisimetrica e das Fracoes Generalizadas de 
    # Autoexpressao de A e B.
    
    # Inicializa motor grafico de plotagem
    plotly()

    # Valores singulares e ordem dos componentes pos ordenamento
    alphas, betas, gsv_sorted_idx = alphas_betas(D1, D2, true)

    # Calcula metricas
    theta = antisymmetric_angular_distance(alphas./betas)
    P1,P2 = generalized_fractions_eigenexpression(alphas, betas)

    # Define os marcadores do eixo X da Distancia Angular Antisimetrica
    ticks = [-pi/4, -pi/8, -3*pi/32, -pi/16, -pi/32, 
                0, pi/32, pi/16, 3*pi/32, pi/8, pi/4]
    tickslabels = ["-pi/4", "-pi/8", "-3pi/32", "-pi/16", "pi/32", 
                  "0", "pi/32", "pi/16", "3pi/32", "pi/8", "pi/4"]

    # Plot da Distancia Angular Antisimetrica 
    p1 = bar(reverse(theta),
      orientation='h',
      title = "Significância dos Componentes de A Relativos a B",
      xlabel = "Distância Angular",
      ylabel = "Componentes",
      xlims = (minimum(ticks)-0.01, maximum(ticks)+0.01),
      xticks = (ticks, tickslabels),
      #ylims = (0,length(theta)+1),
      #yticks = (1:length(theta), [string(i) for i in 1:length(theta)]),
      #ylims = (0,length(gsv_sorted_idx)+1),
      #yticks = (gsv_sorted_idx, [string(i) for i in gsv_sorted_idx]),
      xrotation = 60,
      leg = false,
      titlefontsize = 7,
      title_location = :center,
      tickfontsize = 6,
      legendfontsize = 6,
      guidefontsize = 6,
      legendtitlefontsize = 6
    )

    # Plot da Fracao Generalizada de Autoexpressao de A
    p2 = bar(reverse(P1),
      orientation='h',
      title = "Frações de Autoexpressão para amostras de A",
      xlabel = "Fração de Autoexpressão",
      #ylabel = "Componentes",
      xlims = (+0.01),
      #ylims = (0,length(P1)+1),
      #yticks = (1:length(P1), [string(i) for i in 1:length(P1)]),
      #ylims = (0,length(gsv_sorted_idx)+1),
      #yticks = (gsv_sorted_idx, [string(i) for i in gsv_sorted_idx]),
      xrotation = 60,
      leg = false,
      titlefontsize = 7,
      title_location = :center,
      tickfontsize = 6,
      legendfontsize = 6,
      guidefontsize = 6,
      legendtitlefontsize = 6
    )

    # Plot da Fracao Generalizada de Autoexpressao de B
    p3 = bar(reverse(P2),
      orientation='h',
      title = "Frações de Autoexpressão para amostras de B",
      xlabel = "Fração de Autoexpressão",
      ylabel = "Componentes",
      xlims = (+0.01),
      #ylims = (0,length(P2)+1),
      #yticks = (1:length(P2), [string(i) for i in 1:length(P2)]),
      #ylims = (0,length(gsv_sorted_idx)+1),
      #yticks = (gsv_sorted_idx, [string(i) for i in gsv_sorted_idx]),
      xrotation = 60,
      leg = false,
      titlefontsize = 7,
      title_location = :center,
      tickfontsize = 6,
      legendfontsize = 6,
      guidefontsize = 6,
      legendtitlefontsize = 6
    )

    # Define leiaute de visualizacao dos plots 
    l1 = @layout[a ; b c]

    # Plota metricas
    plot(p1, dpi = 600)
    png(string(output_path, "angular_distance.png"))

    plot(p3, dpi = 600)
    png(string(output_path, "fractions_B.png"))

    plot(p2, dpi = 600)
    png(string(output_path, "fractions_A.png"))

    plot(p1, p3, p2, layout=l1, dpi = 600)
    png(string(output_path, "measures.png"))
end

# ------------------------------------------------------------------------------
# Funcao para localizar exemplos com alta e baixa variancias por componente
# ------------------------------------------------------------------------------

# TODO - Refatorar
function examples_relation_feature_space(A,B)
    #
  
    # Decomposicao GSVD
    U,V,Q,D1,D2,R = svd(A,B)
  
    println("U: ",U)
    println("V: ",V)
  
    # Valores singulares e generalizados
    gsv = Array(diag(D1./D2))
    #svA = Array(diag(D1))
    #svB = Array(diag(D2))
  
    # Quantidade de valores singulares generalizados
    n = length(gsv)
  
    # Fracoes generalizadas de autoexpressao
    #P1,P2 = generalized_fractions_eigenexpression(A,B)
  
    # Inicializa vetor de exemplos por componente
    examples=[]
  
    # Exemplos com maior e menor expressividade, componente a componente
    for i in 1:n
  
      # Indices dos maiores e menores fenomenos por componente e dataset
      min_B = argmin(V[:,i])
      max_B = argmax(V[:,i])
      min_A = argmin(U[:,i])
      max_A = argmax(U[:,i])
  
      # Calcula vetor de exemplos por componente
      push!(examples,[B[max_B,i],B[min_B,i],A[min_A,i],A[max_A,i]])
  
    end
  
    return examples
end
  
# ---------------------------------------------------------------------------
# Funcao para visualizacoes de exemplos nos componentes do GSVD
# ---------------------------------------------------------------------------

# TODO - Refatorar
function show_examples_GSVD(A,B)
    #
  
    plotly()
  
    ticks = [-pi/4, -pi/8, -3*pi/32, -pi/16, -pi/32, 0, pi/32, pi/16, 3*pi/32, pi/8, pi/4]
    tickslabels = ["-pi/4", "-pi/8", "-3pi/32", "-pi/16", "pi/32", "0", "pi/32", "pi/16", "3pi/32", "pi/8", "pi/4"]
  
    # Decomposicao GSVD
    U,V,Q,D1,D2,R = svd(A,B)
    # Valores singulares e generalizados
    alphas, betas = alphas_betas(D1,D2)
    gsv = alphas./betas

    theta = antisymmetric_angular_distance(gsv)
    P1,P2 = generalized_fractions_eigenexpression(alphas,betas)
    gsv_sorted_idx = sortperm(gsv, rev=true)
    examples = examples_relation_feature_space(A,B)

    p1 = bar(theta,
      orientation='h',
      title = "Significância dos Componentes de A Relativos a B",
      xlabel = "Distância Angular",
      ylabel = "Componentes",
      xlims = (minimum(ticks)-0.01, maximum(ticks)+0.01),
      xticks = (ticks, tickslabels),
      #ylims = (0,length(theta)+1),
      #yticks = (1:length(theta), [string(i) for i in 1:length(theta)]),
      ylims = (0,length(gsv_sorted_idx)+1),
      yticks = (gsv_sorted_idx, [string(i) for i in gsv_sorted_idx]),
      xrotation = 60,
      leg = false,
      titlefontsize = 7,
      title_location = :center,
      tickfontsize = 6,
      legendfontsize = 6,
      guidefontsize = 6,
      legendtitlefontsize = 6
    )
    for i in 1:length(examples)
      annotate!(
       [(minimum(ticks)*3/4, i,Plots.text(string(
                round(examples[i][1],digits=4)), :center)),
        (minimum(ticks)*1/4, i,Plots.text(string(
                round(examples[i][2],digits=4)), :center)),
        (maximum(ticks)*1/4, i,Plots.text(string(
                round(examples[i][3],digits=4)), :center)),
        (maximum(ticks)*3/4, i,Plots.text(string(
                round(examples[i][4],digits=4)), :center))]
         )
    end
    p1
end