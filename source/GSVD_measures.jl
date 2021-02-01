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

# ---------------------------------------------------------------------------
# Funcao para escolher os melhores valores singulares de A e B
# ---------------------------------------------------------------------------
# TODO - Elaborar modos 2 e 3 da função
function best_sv_idx(mode, alpha, betas, theta, P1, P2, theta_lim = pi/4)
    # Dado um modo de atuação, os valores singulares do GSVD de A e B, suas 
    # distancias angulares, respectivas autoexpressoes e um limite para 
    # distancia angular, fornece uma lista de indices ordenados dos 
    # "melhores" valores singulares para escolher pinceis. Pinceis são as 
    # linhas da matriz compartilhada da decomposição do GSVD.
    # Os modos da funcao direcionam a definicao de melhor:
    #
    # Modo 1:  Busca componentes com maior autoexpressao e distancia angular
    # aceitavel;
    #
    # Modo 2:
    #
    # Modo 3: 

    if mode == 1
   
        # Lista de distâncias angulares que indiquem similaridade entre os 
        # datasets, compreendidos tal que -theta_lim <= x <= theta_lim.
        sv_ad_list = filter(x -> ((x >= -theta_lim) && (x <= theta_lim)), theta)
        #println("sv_ad_list: ", sv_ad_list) # apenas para debug

        # Lista de índices das distâncias angulares filtradas
        sv_ad_idx_list = indexin(sv_ad_list, theta)
        #println("sv_ad_idx_list: ", sv_ad_idx_list) # apenas para debug

        # Lista de fracoes de autoexpressao, correspondentes as distancias 
        # angulares filtradas. 
        sv_ae_list = getindex(P1, sv_ad_idx_list)
        #println("sv_ae_list: ", sv_ae_list) # apenas para debug

        # Ordena lista de fracoes de autoexpressao
        sv_ae_list = sort(sv_ae_list, rev = true)
        #println("sv_ae_list ordened: ", sv_ae_list) # apenas para debug

        # Lista de índices das fracoes de autoexpressao ordenadas
        sv_ae_idx_list = indexin(sv_ae_list, P1)
        #println("sv_ae_idx_list: ", sv_ae_idx_list) # apenas para debug

        # FIXME - As frações de autoexpressao têm que considerar P2. Provavelmente buscar algum P1/P2 próximo a 1.
        
        # FIXME - Função indexin() retorna o índice do primeiro valor que faz match, o que pode ocasionar um range externo ao filtrado.
        # Ex.:  Sejam   thetas          = [-20, -10, 0, 10, 20, 30]
        #               theta_lim       = 10
        #               alphas          = [0.2, 0.2, 0.05, 0.4, 0.05, 0.3]
        #
        # É retornado   sv_ad_list      = [-10, 0, 10]
        #               sv_ad_idx_list  = [2, 3, 4]
        #               sv_ae_list      = [0.4, 0.2, 0.05]
        #               sv_ae_idx_list  = [4, 1, 3]
        #
        # Analisando as variáveis sv_ad_idx_list e sv_ae_idx_list, percebe-se que o índice 1 não pertence a ambos índices.

    elseif mode == 2
        #

    elseif mode == 3
        #

    end

    # retorna lista de indices
    return sv_ae_idx_list 
end

# TODO - Criar testes das funcoes