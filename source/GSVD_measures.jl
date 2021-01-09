# ---------------------------------------------------------------------------
# Bibliotecas
# ---------------------------------------------------------------------------

using Random # necessario para funcoes sort


# ---------------------------------------------------------------------------
# Funcao para Valores Singulares de A e B
# ---------------------------------------------------------------------------

function alphas_betas(D1, D2, idx = false)
    # Dadas as matrizes D1 e D2 da fatoracao GSVD(A,B), retorna os valores
    # singulares (alphas e betas) das matrizes A e B, respectivamente.
    # Opcionalmente, e possivel retornar a ordem dos componentes pos 
    # ordenamento de alphas e betas.
    

    # Como os valores singulares na matriz D2 nem sempre estão localizados na 
    # primeira diagonal da respectiva, é preciso extrai-la na posicao correta.
    
    # Localiza onde inicia a diagonal em D2
    k = argmax(D2[1, :]) - 1
    
    # Extrai os valores singulares das respectivas diagonais
    alphas = Array(diag(D1))
    betas  = Array(diag(D2, k))


    # Como os valores singulares na matriz D2 nem sempre estão localizados na 
    # primeira diagonal da respectiva, é preciso incluir os valores nulos, 
    # garantindo a mesma quandidade de alfas e betas, que é igual ao tamanho 
    # da dimensao compartilhada das matrizes A e B.

    # Calcula vetores de valores nulos. As quantidades se baseiam na diferenca
    # entre o tamanho da dimensao compartilhada de A e B e o quantidade de 
    # valores singulares presentes nas diagonais de D1 e D2.
    zeros_alpha = zeros(size(D1)[2] - length(alphas))
    zeros_beta  = zeros(size(D2)[2] - length(betas ))
    
    # Atualiza valores singulares com os valores nulos
    alphas = vcat(alphas    , zeros_alpha)
    betas  = vcat(zeros_beta, betas      )

    # Calcula os indices para ordenamento dos valores singulares, baseando-se 
    # no ordenamento dos valores singulares generalizados
    gsv_sorted_idx = sortperm(alphas ./ betas, rev = true)
  
    # Ordenamento dos valores singulares
    alphas = getindex(alphas, gsv_sorted_idx)
    betas  = getindex(betas , gsv_sorted_idx)
    
    # Verifica solicitacao do retorno da ordem dos componentes, pos
    # ordenamento de alphas e betas.
    if idx == true
        return alphas, betas, gsv_sorted_idx
    else
        return alphas, betas
    end
end


# ---------------------------------------------------------------------------
# Funcao para Distancia Angular Antisimetrica
# ---------------------------------------------------------------------------

function antisymmetric_angular_distance(gsv, idx = false)
    # Dados os valores singulares generalizados do GSVD(A,B), retorna a 
    # distancia angular antisimetrica (theta).
    # Opcionalmente, e possivel retornar a ordem dos componentes pos 
    # ordenamento dos valores singulares.

    # Calcula os indices para ordenamento dos valores singulares, baseando-se 
    # no ordenamento dos valores singulares generalizados
    gsv_sorted_idx = sortperm(gsv, rev = true)

    # Ordenamento dos valores singulares generalizados
    gsv = getindex(gsv, gsv_sorted_idx)
  
    # Inicializa vetor de distancias angulares antisimetricas
    theta = Float64[]
  
    # Calcula o vetor de distancias angulares antisimetricas
    for i in 1:length(gsv)
      push!(theta, atan(gsv[i]) - (pi / 4))
    end
    
    # Verifica solicitacao do retorno da ordem dos componentes, pos
    # ordenamento de alphas e betas.
    if idx == true
        return theta, gsv_sorted_idx
    else
        return theta
    end
end


# ---------------------------------------------------------------------------
# Funcao para Fracoes Generalizadas de Autoexpressao
# ---------------------------------------------------------------------------

function generalized_fractions_eigenexpression(alphas, betas, idx = false)
    # Dados os valores singulares (alphas e betas) das matrizes A e B do 
    # GSVD(A,B), retorna as fracoes generalizadas de autoexpressao (P1 e P2).
    # Opcionalmente, e possivel retornar a ordem dos componentes pos 
    # ordenamento dos valores singulares.
  
    # Calcula os indices para ordenamento dos valores singulares, baseando-se 
    # no ordenamento dos valores singulares generalizados
    gsv_sorted_idx = sortperm(alphas ./ betas, rev = true)
  
    # Ordenamento dos valores singulares
    alphas = getindex(alphas, gsv_sorted_idx)
    betas  = getindex(betas , gsv_sorted_idx)

    # Calcula as fracoes generalizadas de autoexpressao 
    P1 = (alphas.^ 2) ./ (sum(alphas.^ 2))
    P2 = (betas .^ 2) ./ (sum(betas .^ 2))
    
    # Verifica solicitacao do retorno da ordem dos componentes, pos
    # ordenamento de alphas e betas.
    if idx == true
        return P1, P2, gsv_sorted_idx
    else
        return P1, P2
    end
end

# TODO - Criar testes das funcoes