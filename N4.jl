
using Plots
using TaylorSeries 

#Esa funcion calculará los primeros coeficintes                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       #Definiremos una función para los primeros datos
function PrimerosD(a_1,epsilon)
    a_0=b_0=beta_0=0.0
    alpha_0=1.0
    a=Float64[]  #Especificamos el tipo de bichos que tendrá la lista
    b=Float64[]
    beta=Float64[]
    alpha=Float64[]
    
    push!(a,a_0)  #Agregamos el primer elemento a cada lista
    push!(b,b_0)
    
    push!(alpha,alpha_0)
    push!(beta,beta_0)

    push!(a,a_1)
  # usando el parametro epsilon calculamos los siguientes 2 elementos
    l=1.0+epsilon/2.+((epsilon*(4+epsilon))^.5)/2.
    b_1 = a_1*(.5+(4.*epsilon+epsilon^2)^.5/(2*epsilon))
    push!(b,b_1)
    
    #  Actualizamos alpha y beta
    alpha_1=-beta_0*b_1
    push!(alpha,alpha_1)
    beta_1=alpha_0*b_1
    push!(beta,beta_1)
    return a,b,alpha,beta,l

end


# nn equivale a la "n" en el artículo.
# Alpha y betha son las listas que guardan a las mismas variables.
function sigmaa(nn, alpha, beta)
    SumSigmaA = 0.0 
    for k = 0:nn-2
        SigmaA =(k+1)*alpha[nn-k]*beta[k+2]
        SumSigmaA = SigmaA+SumSigmaA
    end
    SumSigmaA
end



function sumaal(nn,beta,b)
    SumA = 0.0
    for k = 0:nn-1
        SA = (k+1)*beta[nn-k]*b[k+2]
        SumA = SA + SumA
    end
    SumA
end

#Definimos una función que va a calcular todos los términos de lo polinomios,
# solo recibe los parametros iniciales y regresa dos listas y dos polinomios, 
# los polinomios que regresa tienen los coeficientes de la lista
function Polinomio(o,a_1,p)
    a,b,alpha,beta,l=PrimerosD(a_1,p)
    for n = 3:o
        SumSigmaA = 0.0
        SumSigmaB = 0.0
        SumA = 0.0
        SumB = 0.0
        nn=n-1
    
    SumSigmaA = sigmaa(nn, alpha, beta)
    a_n = ((-p*(1-l^(n+1)))/((n+1)*((1-l^(n+1))*(1+p-l^(n+1))-p)))*SumSigmaA
    push!(a,a_n)
   
    SumSigmaB = sigmaa(nn, alpha, b)
    b_n = ((p*(l^(n+1)))/((n+1)*((1-l^(n+1))*(1+p-l^(n+1))-p)))*SumSigmaB
    push!(b,b_n)
    
    SumA=sumaal(nn,beta,b)
    alpha_n = ((-1)/(n+1))*SumA
    push!(alpha,alpha_n)
    
    SumB=sumaal(nn,alpha,b)
    beta_n = (1/(n+1))*SumB
    push!(beta,beta_n)
    
    
    end
    A=Taylor1(a)
    B=Taylor1(b)
    return a,b,A,B
end
    


#con esta funcion ya tenemos dos listas y dos polinomios
a,b,A,B=Polinomio(10,0.2,0.3)   
print(a,A)



#aquí lo que hacemos es evaluar el polinomio en un intervalo de tiempo
ValX=Float64[]
ValY=Float64[]

for t = 0:0.125:17.
    x = mod2pi( evaluate(A,t))
    y = mod2pi(evaluate(B,t))
    push!(ValX,x)
    push!(ValY,y)
end

plot(ValX,ValY)

#definimos la función del mapeo, la cual recibe dos polinomios y les aplica el mapeo
# como resultado regresa dos polinomios.
function Mapeo(X,Y,ϵ)
    x=X+ϵ*sin(Y)
    y=X+Y+ϵ*sin(Y)
    return x,y
end

#definimos la funcion para el polinomio, la cual recibe dos listas y les aplica la 
# composicion de funciones de manera directa, al final usa las listas para devolver un
# polinomio.
function Pol(X,Y,ϵ)
    x=Float64[]
    y=Float64[]
    λ=1.0+ϵ/2.+((ϵ*(4+ϵ))^.5)/2.
    for i=1:length(X)
        push!(x,X[i]*λ^(i-1))
        push!(y,Y[i]*λ^(i-1))
    
    end
    x=Taylor1(x)
    y=Taylor1(y)
    return x,y
end
    
        

#ya que tenemos las dos composiciones cada una con dos polinomios podemos calcular el error,
# aquí la función Err recibe como entrada los 4 polinomios, el tiempo y el paso del tiempo,
# regresa una gráfica con los valores del error definido como ya s explico, se usa la norma
# infinito para ello. 
function Err(P1,P2,K1,K2,T,paso)
    Val=Float64[]
    Tiem=Float64[]
    E1=P1-K1
    E2=P2-K2
    for t = 0:paso:T
        x = evaluate(E1,t)
        y = evaluate(E2,t)
        E=[x,y]
        norma = norm(E,Inf)
    
        push!(Val,norma)
        push!(Tiem,t)
    
    end
    plot!(Tiem,Val)
end

#Esta es una función que calcula todos los pasos anteriores, dado un polinomio
# calcula las dos composicioes y luego los evalúa para encntrar el error.
function ErrorE(a,b,A,B,ϵ,T,paso)
    
    Xn,Yn=Mapeo(A,B,ϵ)
    Xnn,Ynn=Pol(a,b,ϵ)
    Err(Xn,Yn,Xnn,Ynn,T,paso)
end

A,B,a,b=Polinomio(10,0.2,0.3)   
X,Y,x,y=Polinomio(12,0.2,0.3)
W,Z,w,z=Polinomio(14,0.2,0.3)
Q,R,q,r=Polinomio(16,0.2,0.3)


ErrorE(A,B,a,b,0.3,12,0.025)
ErrorE(X,Y,x,y,0.3,12,0.025)
ErrorE(W,Z,w,z,0.3,12,0.025)
ErrorE(Q,R,q,r,0.3,12,0.025)
    


