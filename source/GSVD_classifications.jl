# ---------------------------------------------------------------------------
# Bibliotecas
# ---------------------------------------------------------------------------

using LinearAlgebra


# ---------------------------------------------------------------------------
# Funcao para criar matrizes ortogonais
# ---------------------------------------------------------------------------

function orth(matrix)
    #Dada uma matriz, retorna uma base orthogonal de mesma dimensao
    
    # Decomposicao QR da matriz, onde Q é uma base ortogonal
    Q, = qr(matrix)
    
    return Q
end


# ---------------------------------------------------------------------------
# Funcao para criar relacoes
# ---------------------------------------------------------------------------

# FIXME - Criação de matrizes não respeita posição e dimensões das decomposições de matrizes gordas
function make_relation(m,p,n,svA,svB)
    # Dadas as dimensoes das matrizes A (mxn) e B (pxn), alem dos seus 
    # respectivos vetores de valores singulares, retorna uma relacao baseada 
    # na decomposicao GSVD.
    
    # Matrizes ortogonais aleatorias
    U=orth(rand(m, n))
    V=orth(rand(p, n))
    
    # Matrizes diagonais de valores singulares
    D1=diagm(svA)
    D2=diagm(svB)
    
    # Matriz aleatoria
    RQ=rand(n, n)
    
    # Reconstrucao das matrizes baseado na decomposicao GSVD
    A=U*D1*RQ'
    B=V*D2*RQ'
    
    # Valores singulares de A e B
    println("svA: ",svA)
    println("svB: ",svB)
    
    return A,B
end


# ---------------------------------------------------------------------------
# Funcao para criar relacoes iguais
# ---------------------------------------------------------------------------

# FIXME - Criação de matrizes não respeita posição e dimensões das decomposições de matrizes gordas
function make_equ_relation(m, p, n)
    #D adas as dimensoes das matrizes A (mxn) e B (pxn), retorna 2 matrizes com uma relacao igual
    
    # Valores singulares aleatorios de A
    singular_values_A = rand(n)
    
    # Valores singulares de B iguais aos de A
    singular_values_B = singular_values_A
    
    return make_relation(m, p, n, singular_values_A, singular_values_B)
end


# ---------------------------------------------------------------------------
# Funcao para criar relacoes diferentes
# ---------------------------------------------------------------------------

# FIXME - Criação de matrizes não respeita posição e dimensões das decomposições de matrizes gordas
function make_dif_relation(m, p, n, lim_dif=50.02)
    # Dadas as dimensoes das matrizes A (mxn) e B (pxn), retorna 2 matrizes com uma relacao de diferenca, sendo possivel definir um limite para a diferenca
    
    # Inicializa vetores de valores singulares
    singular_values_A=Float64[]
    singular_values_B=Float64[]
    
    # Partes dos valores singulares
    div = rand(1:n-1)
    
    # Parte dos valores singulares de A maiores que os de B
    for i in 1:(div)
      
      # Valor singular de B aleatorio
      svB=rand()
      
      # A razão entre os valores deve ser >= limite da diferenca
      svA=svB*(lim_dif)
      
      # Atualiza vetores de valores singulares
      push!(singular_values_A,svA) 
      push!(singular_values_B,svB)
    end
    
    # Parte dos valores singulares de A menores que os de B
    for i in (div+1):n
      
      # Valor singular de B aleatorio
      svB=rand()
      
      # A razão entre os valores deve ser <= 1/limite da diferenca
      svA=svB*(1/lim_dif)
      
      # Atualiza vetores de valores singulares
      push!(singular_values_A,svA)
      push!(singular_values_B,svB)
    end
    
    return make_relation(m, p, n, singular_values_A, singular_values_B)
end


# ---------------------------------------------------------------------------
# Funcao para criar relacoes de contingencia
# ---------------------------------------------------------------------------

# FIXME - Criação de matrizes não respeita posição e dimensões das decomposições de matrizes gordas
function make_con_relation(m, p, n, kind, lim_equ=1.02, lim_dif=50.02)
    # Dadas as dimensoes das matrizes A (mxn) e B (pxn) e o tipo de contingencia (contido/contem), retorna 2 matrizes com uma relacao de contingencia, sendo possivel definir um limite para a igualdade e diferenca
   
    # Inicializa vetores de valores singulares
    singular_values_A=Float64[]
    singular_values_B=Float64[]
    
    # Partes dos valores singulares
    div = rand(1:n-1)
    
    # Tipo A contem B
    if kind == "contains"
      
      # Parte dos valores singulares de A maiores que os de B
      for i in 1:(div)
        
        # Valor singular de B aleatorio
        svB=rand()
        
        # A razão entre os valores deve ser >= limite da diferenca
        svA=svB*(lim_dif)
        
        # Atualiza vetores de valores singulares
        push!(singular_values_A,svA)
        push!(singular_values_B,svB)
      end
      
      # Parte dos valores singulares de A iguais aos de B
      for i in (div+1):n
        
        # Valor singular de B aleatorio
        svB=rand()
        
        # A razão entre os valores deve estar compreendido no limite da igualdade
        svA=svB
        
        # Atualiza vetores de valores singulares
        push!(singular_values_A,svA)
        push!(singular_values_B,svB)
      end
    end
    
    # Tipo A esta contido em B
    if kind == "contained"
      
      # Parte dos valores singulares de A menores que os de B
      for i in 1:(div)
        
        # Valor singular de B aleatorio
        svB=rand()
        
        # A razão entre os valores deve ser <= 1/limite da diferenca
        svA=svB*(1/lim_dif)
        
        # Atualiza vetores de valores singulares
        push!(singular_values_A,svA)
        push!(singular_values_B,svB)
      end
      
      # Parte dos valores singulares de A iguais aos de B
      for i in (div+1):n
        
        # Valor singular de B aleatorio
        svB=rand()
        
        # A razão entre os valores deve estar compreendido no limite da igualdade
        svA=svB
        
        # Atualiza vetores de valores singulares
        push!(singular_values_A,svA)
        push!(singular_values_B,svB)
      end
    end
    
    return make_relation(m, p, n, singular_values_A, singular_values_B)
end


# ---------------------------------------------------------------------------
# Funcao para criar relacoes de intersecao
# ---------------------------------------------------------------------------

# FIXME - Criação de matrizes não respeita posição e dimensões das decomposições de matrizes gordas
function make_int_relation(m, p, n, lim_equ=1.02, lim_dif=50.02)
    # Dadas as dimensoes das matrizes A (mxn) e B (pxn), retorna 2 matrizes com uma relacao de intersecao, sendo possivel definir um limite para a igualdade e diferenca
    
    # Inicializa vetores de valores singulares
    singular_values_A=Float64[]
    singular_values_B=Float64[]
    
    # Partes dos valores singulares
    div1 = rand(       1:(n-2))
    div2 = rand((div1+1):(n-1))
    
    # Parte dos valores singulares de A maiores que os de B
    for i in 1:(div1)
      
      # Valor singular de B aleatorio
      svB=rand()
      
      # A razão entre os valores deve ser >= limite da diferenca
      svA=svB*(lim_dif)
      
      # Atualiza vetores de valores singulares
      push!(singular_values_A,svA)
      push!(singular_values_B,svB)
    end
    
    # Parte dos valores singulares de A iguais aos de B
    for i in (div1+1):(div2)
      
      # Valor singular de B aleatorio
      svB=rand()
      
      # A razão entre os valores deve estar compreendido no limite da igualdade
      svA=svB
      
      # Atualiza vetores de valores singulares
      push!(singular_values_A,svA)
      push!(singular_values_B,svB)
    end
    
    # Parte dos valores singulares de A menores que os de B
    for i in (div2+1):n
      
      # Valor singular de B aleatorio
      svB=rand()
      
      # A razão entre os valores deve ser <= 1/limite da diferenca
      svA=svB*(1/lim_dif)
      
      # Atualiza vetores de valores singulares
      push!(singular_values_A,svA)
      push!(singular_values_B,svB)
    end
    
    return make_relation(m, p, n, singular_values_A, singular_values_B)
end


# ---------------------------------------------------------------------------
# Funcao para classificar relacoes
# ---------------------------------------------------------------------------

# FIXME - Criação de matrizes não respeita posição e dimensões das decomposições de matrizes gordas
function classify_relation(A, B, lim_equ=1.01, lim_dif=50.01)
    # Dadas as matrizes A e B de uma relacao, retorna a classificacao desta relacao(igualdade,diferenca,contingecia,intersecao), sendo possivel definir um limite para a igualdade e diferenca
    
    # Decomposicao GSVD
    U,V,Q,D1,D2,R = svd(A,B)
    
    # Valores singulares e generalizados
    gsv = Array(diag(D1./D2))
    
    # Quantidade de valores singulares generalizados
    n = length(gsv)
    
    # Classificacao dos tipos dos valores singulares generalizados
    #
    #   high_B_gsv  |  low_B_gsv  |  equal_gsv  |  low_A_gsv  |  high_A_gsv
    #               |             |             |             |             
    # 0        1/lim_dif     1/lim_equ   1   lim_equ       lim_dif          +∞
    # <_____________|_____________|______|______|_____________|_____________>
    
    # Contadores de tipos
    equal_gsv =0
    low_A_gsv =0
    low_B_gsv =0
    high_A_gsv=0
    high_B_gsv=0
    
    # Contando os valores singulares generalizados por tipo
    for i in 1:n
      
      # Tipo "igual"
      # Compreendido entre o limite inferior e superior de igualdade
      ((gsv[i] >= 1/lim_equ) && (gsv[i] <=   lim_equ)) && (equal_gsv  += 1)
      
      # Tipo "fraco para A"
      # Compreendido entre o limite superior de igualdade e superior da diferenca
      ((gsv[i] >    lim_equ) && (gsv[i] <    lim_dif)) && (low_A_gsv  += 1)
      
      # Tipo "fraco para B"
      # Compreendido entre o limite inferior da diferenca e inferior da igualdade
      ((gsv[i] >  1/lim_dif) && (gsv[i] <  1/lim_equ)) && (low_B_gsv  += 1)
      
      # Tipo "forte para A"
      # Maior que o limite superior da diferenca
       (gsv[i] >=   lim_dif)                           && (high_A_gsv += 1)
      
      # Tipo "forte para B"
      # Maior que o limite inferior da diferenca
       (gsv[i] <= 1/lim_dif)                           && (high_B_gsv += 1)
    
    end
    
    # Classificacao de relacoes de acordo com os tipos de valores singulares
    
    # Necessário apenas tipos "igual"
    if (equal_gsv > 0) && (equal_gsv == n)
      return classif = "As matrizes A e B são iguais!"
    end
    
    # Necessário apenas tipos "forte para A" e "forte para B"
    if (high_A_gsv > 0) && (high_B_gsv > 0) && (high_A_gsv + high_B_gsv == n)
      return classif = "As matrizes A e B são diferentes!"
    end
    
    # Necessários tipos "forte para A" ou "fraco para A", aceitando tipo "igual"
    if ((high_A_gsv > 0) || (low_A_gsv > 0)) && (high_A_gsv + low_A_gsv + equal_gsv == n)
      return classif = "A matriz A contém a matriz B!"
    end
    
    # Necessários tipos "forte para B" ou "fraco para B", aceitando tipo "igual"
    if ((high_B_gsv > 0) || (low_B_gsv > 0)) && (high_B_gsv + low_B_gsv + equal_gsv == n)
      return classif = "A matriz A está contida na matriz B!"
    end
    
    # Necessários tipos "forte para A" ou "fraco para A", e "forte para B" ou "fraco para B" e "igual"
    if ((high_A_gsv > 0) || (low_A_gsv > 0)) && ((high_B_gsv > 0) || (low_B_gsv > 0)) && (equal_gsv > 0)
      return classif = "As matrizes A e B possuem uma interseção!"
    end
    
    # Necessários tipo "fraco para A" e "fraco para B"
    if (low_A_gsv > 0) && (low_B_gsv > 0) && (low_A_gsv + low_B_gsv == n)
      return classif = "Não sabemos definir essa classificação"
    end
    
    return classif
end

# ------------------------------------------------------------------------------
# Funcao para classificar relacoes com fracoes de autoexpressao
# ------------------------------------------------------------------------------

# TODO - Finalizar funcao classificadora com autoexpressoes
function classify_relation_v2(A, B, lim_equ=1.01, lim_dif=50.01, lim_aex=0.02)
    # Dadas as matrizes A e B de uma relacao, retorna a classificacao desta relacao(igualdade,diferenca,contingecia,intersecao), sendo possivel definir um limite para a igualdade, diferenca e autoexpressao
  
    # Decomposicao GSVD
    U,V,Q,D1,D2,R = svd(A,B)
  
    # Valores singulares e generalizados
    gsv = Array(diag(D1./D2))
  
    # Quantidade de valores singulares generalizados
    n = length(gsv)
  
    # Fracoes generalizadas de autoexpressao
    P1,P2 = generalized_fractions_eigenexpression(A,B)
  
    # Valores singulares e generalizados
    #D1 = Array(diag(D1))
    #D2 = Array(diag(D2))
    #gsv = Array((D1./D2))
  
    # Classificacao dos tipos das fracoes generalizadas de autoexpressao
    #
    #         low_A_gfe/low_B_gfe        |        high_A_gfe/high_B_gfe
    #                                    |
    # 0                               lim_aex                               1
    # <__________________________________|__________________________________>
  
    # Classificacao dos tipos dos valores singulares generalizados
    #
    #   high_B_gsv  |  low_B_gsv  |  equal_gsv  |  low_A_gsv  |  high_A_gsv
    #               |             |             |             |
    # 0        1/lim_dif     1/lim_equ   1   lim_equ       lim_dif          +∞
    # <_____________|_____________|______|______|_____________|_____________>
  
    # Contadores de tipos
    equal_gsv =0
    low_A_gsv =0
    low_B_gsv =0
    high_A_gsv=0
    high_B_gsv=0
    low_A_gfe =0
    low_B_gfe =0
    high_A_gfe=0
    high_B_gfe=0
  
      # Contando as fracoes generalizadas de autoexpressao e os valores singulares generalizados por tipo
    for i in 1:n
  
      # Fracoes generalizadas de autoexpressao
  
      # Tipo "gfe baixo para A"
      # Menor ou igual ao limite da autoexpressao
      (P1[i] <= lim_aex)                               && (low_A_gfe  += 1)
  
      # Tipo "gfe baixo para B"
      # Menor ou igual ao limite da autoexpressao
      (P2[i] <= lim_aex)                               && (low_B_gfe  += 1)
  
      # Tipo "gfe alto para A"
      # Maior que o limite da autoexpressao
      (P1[i] >  lim_aex)                               && (high_A_gfe += 1)
  
      # Tipo "gfe alto para B"
      # Maior que o limite da autoexpressao
      (P2[i] >  lim_aex)                               && (high_B_gfe += 1)
  
      # Valores singulares generalizados
  
      # Tipo "gsv igual"
      # Compreendido entre o limite inferior e superior de igualdade
      ((gsv[i] >= 1/lim_equ) && (gsv[i] <=  lim_equ) && (P1[i] >  lim_aex) && (P2[i] >  lim_aex)) && (equal_gsv  += 1)
  
      # Tipo "gsv fraco para A"
      # Compreendido entre o limite superior de igualdade e superior da diferenca
      ((gsv[i] >    lim_equ) && (gsv[i] <   lim_dif) && (P1[i] >  lim_aex) && (P2[i] >  lim_aex)) && (low_A_gsv  += 1)
  
      # Tipo "gsv fraco para B"
      # Compreendido entre o limite inferior da diferenca e inferior da igualdade
      ((gsv[i] >  1/lim_dif) && (gsv[i] < 1/lim_equ) && (P1[i] >  lim_aex) && (P2[i] >  lim_aex)) && (low_B_gsv  += 1)
  
      # Tipo "gsv forte para A"
      # Maior que o limite superior da diferenca
      ((gsv[i] >=   lim_dif) && (P1[i] >  lim_aex) && (P2[i] >  lim_aex)) && (high_A_gsv += 1)
  
      # Tipo "gsv forte para B"
      # Maior que o limite inferior da diferenca
      ((gsv[i] <= 1/lim_dif) && (P1[i] >  lim_aex) && (P2[i] >  lim_aex)) && (high_B_gsv += 1)
  
    end
  
    # Classificacao de relacoes de acordo com os tipos de valores singulares e tipos de fracoes generalizadas
  
    classif = "Não sabemos definir essa classificação"
  
    # Necessário apenas tipos "gsv igual"
    if (equal_gsv > 0) && (equal_gsv - low_A_gfe - low_B_gfe == n)
      return classif = "As matrizes A e B são iguais!"
    end
  
    # Necessário apenas tipos "gsv forte para A" e "gsv forte para B"
    if (high_A_gsv > 0) && (high_B_gsv > 0) && (high_A_gsv + high_B_gsv  - low_A_gfe - low_B_gfe == n)
      return classif = "As matrizes A e B são diferentes!"
    end
  
    # Necessários tipos "gsv forte para A" ou "gsv fraco para A", aceitando tipo "gsv igual"
    if ((high_A_gsv > 0) || (low_A_gsv > 0)) && (high_A_gsv + low_A_gsv + equal_gsv - low_A_gfe - low_B_gfe == n)
      return classif = "A matriz A contém a matriz B!"
    end
  
    # Necessários tipos "gsv forte para B" ou "gsv fraco para B", aceitando tipo "gsv igual"
    if ((high_B_gsv > 0) || (low_B_gsv > 0)) && (high_B_gsv + low_B_gsv + equal_gsv - low_A_gfe - low_B_gfe == n)
      return classif = "A matriz A está contida na matriz B!"
    end
  
    # Necessários tipos "gsv forte para A" ou "gsv fraco para A", e "gsv forte para B" ou "gsv fraco para B" e "gsv igual"
    if ((high_A_gsv > 0) || (low_A_gsv > 0)) && ((high_B_gsv > 0) || (low_B_gsv > 0)) && (equal_gsv > 0) && (high_B_gsv + low_B_gsv + equal_gsv + low_A_gsv + high_A_gsv - low_A_gfe - low_B_gfe == n)
      return classif = "As matrizes A e B possuem uma interseção!"
    end
  
    # Necessários tipo "gsv fraco para A" e "gsv fraco para B"
    if (low_A_gsv > 0) && (low_B_gsv > 0) && (low_A_gsv + low_B_gsv - low_A_gfe - low_B_gfe  == n)
      return classif = "Não sabemos definir essa classificação"
    end
  
    return classif
end


# ---------------------------------------------------------------------------
# Testes
# ---------------------------------------------------------------------------

function test_classifications(m,p,n)
  # Dadas as dimensoes das matrizes A (mxn) e B (pxn), retorna testes das 
  # diferentes classificacoes de relacoes possiveis.

  # Relacao de igualdade aleatoria
  println("Criada relação A(",m,"x",n,") igual à B(",p,"x",n,").")
  A, B = make_equ_relation(m,p,n)
  println("Resposta do classificador: ",classify_relation(A,B),"\n")

  # Relacao de diferenca aleatoria
  println("Criada relação A(",m,"x",n,") diferente de B(",p,"x",n,").")
  A, B = make_dif_relation(m,p,n)
  println("Resposta do classificador: ",classify_relation(A,B),"\n")

  # Relacao de contingencia aleatoria, onde A contem B
  println("Criada relação A(",m,"x",n,") contém B(",p,"x",n,").")
  A, B = make_con_relation(m,p,n,"contains")
  println("Resposta do classificador: ",classify_relation(A,B),"\n")

  # Relacao de contingencia aleatoria, onde A esta contido em B
  println("Criada relação A(",m,"x",n,") está contido em B(",p,"x",n,").")
  A, B = make_con_relation(m,p,n,"contained")
  println("Resposta do classificador: ",classify_relation(A,B),"\n")

  # Relacao de intersecao aleatoria
  println("Criada relação A(",m,"x",n,") possui interseção com B(",p,"x",n,").")
  A, B = make_int_relation(m,p,n)
  println("Resposta do classificador: ",classify_relation(A,B),"\n")

end