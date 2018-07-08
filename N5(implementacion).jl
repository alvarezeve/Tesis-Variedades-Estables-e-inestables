																																																																													#usaremos el paquete de TaylorSeries
using TaylorSeries
#definimos la función del mapeo 
function StandarMap(x,p,k)
    x_n = x+k*sin(p)
    p_n = x_n + p
    return [x_n,p_n]
end
#Sabemos que en este caso el punto fijo es 
p_0=[0,0]

#--------------ORDEN 1 ---------------------------------------------
#Definimos un polinomio que será la variable muda de las variedades, como un taylor1
t=Taylor1([0,1],5)

#dado que queremos que x, p sean polinomios necesitamos construirlos, 
# a priori no sabemos mas que el primer coeficiente de estos que es el cero pero 
#la idea es ir calculando los términos en cada iteración, por lo que ponemos los coeficientes que buscamos
# como un x, p 
x,p = set_variables(Float64,"x p",order=1)

#Ahora necesitamos los polinomios en específico, que serán los que metamos al mapeo 
tx=evaluate(t,x)
tp=evaluate(t,p)
#Aplicamos el mapeo a estos dos polinomios
PO=StandarMap(tx,tp,0.3)
# le calculamos el jacobiano, lo que nos dará la ecuación matricial cohomologica
JPO=jacobian(PO)
# le calculamos los eigenvalores y eigenvectores. 
eval,eve=eig(JPO)
# dado que usaremos la variedad inestable, entonces usamos el segundo eigenvalor, vector
λ=eval[2]
xn=eve[1,2]
pn=eve[2,2]

vec=[xn,pn]
# estos coeficientes xn pn serán los coeficientes de primer orden en nuestra variedad

#-------------ORDEN 2-------------------------------------------------
#Queremos hacer la misma idea para este orden
t=Taylor1([0,1],3)
x,p = set_variables(Float64,"x p",order=3)
#usamos los valores ya calculados en el orden 1 para x,p de primer orden
tx=Taylor1([0,xn,x],3)
tp=Taylor1([0,pn,p],3)
#aplicamos el mapeo
SO=StandarMap(tx,tp,0.3)
# como nos falta el otro lado de la ecuación cohomológica definimos un vector de lambdas
vλ=[1,λ,λ^2]
#y definimos también un polinomio para las lambdas
xλt=Taylor1([0,xn*λ,x*λ^2],2)
pλt=Taylor1([0,pn*λ,p*λ^2],2)
λvec=[xλt,pλt]
# ahora ya tengo las dos partes de la ecuación y debo igualarlas para resolver.
Ecua=SO-λvec

#---nota, en el primero coeficiente no se anula en la x, sale del orden de 10^-16---
# de esta ecuación necesitamos solo los de segundo orden, así que los extraemos manualmente 
x2=Ecua[1].coeffs[3]
p2=Ecua[2].coeffs[3]
vec2=[x2,p2]
#calculamos ahora el jacobiano de este vector, ya que los coeficientes serán parte de la matriz que queremos
JSO=jacobian(vec2)
#calculamos su determinante para ver si es cero, 
det(JSO)
# como es distinto de cero y la ecuación esta igualada a cero entonces la única solución es la trivial
x2=0.
p2=0.
vec2=[0,0]

#----------ORDEN 3--------------------------------------------------
t=Taylor1([0,1],3)
x,p = set_variables(Float64,"x p",order=3)
tx=Taylor1([0,xn,x2,x],4)
tp=Taylor1([0,pn,p2,P],4)
TO=StandarMap(tx,tp,0.3)
xλt=Taylor1([0,xn*λ,x2*λ^2,x*λ^3],4)
pλt=Taylor1([0,pn*λ,p2*λ^2,p*λ^3],4)
λvec=[xλt,pλt]
Ecua=TO-λvec
x3=Ecua[1].coeffs[4]
p3=Ecua[2].coeffs[4]

#Dado que en este momento comienzan a aparecer términos independientes en x, p necesitamos una forma de extraerlos
# primero estraemos los coeficientes de x, p como antes con el jacobiano. 
JTO=jacobian([x3,p3])
#calculamos el determinante 
det(JTO)
#extraemos ahora los coeficientes independientes de x,p 
a= Ecua[1].coeffs[4].coeffs[1]
b= Ecua[2].coeffs[4].coeffs[1]
vecCoef=[a,b]

#ya que tenemos los coeficientes necesitamos reolver la ecuación JTO[x3,p3]=vecC[x3,p3]
# entonces solo se trata de invertir el jacobiano y multiplicar con el vector del lado izquierdo (¿como quitar que sea pol homogeneo?)

T3=transpose(vecCoef)*(1\JTO)
x3=T3[1]
p3=T3[2]
#Estos dos valores son los coeficientes de x, p de tercer orden en ele polinomo de x,y
#el pequeño problema con estos dos términos es que son tipo TaylorSeries.HomogeneousPolynomial{Float64}, y cuando los queremos usar de 
# de nuevo hay problema al calcular el sin(TaylorSeries.HomogeneousPolynomial{Float64})

#-------------ORDEN 4 -------------------------------------
t=Taylor1([0,1],5)
x,p = set_variables(Float64,"x p",order=5)
tx=Taylor1([0,xn,x2,x3],5)
tp=Taylor1([0,pn,p2,p3],5)
CO=StandarMap(tx,tp,0.3)
xλt=Taylor1([0,xn*λ,x2*λ^2,x3*λ^3,x*λ^4],4)
pλt=Taylor1([0,pn*λ,p2*λ^2,p3*λ^3,p*λ^4],4)
λvec=[xλt,pλt]
Ecua=TO-λvec
x4=Ecua[1].coeffs[5]
p4=Ecua[2].coeffs[5]


JTO=jacobian([x4,p4])







######## .................lo que hoce con Luis
julia> get_coeff(Ecua[1].coeffs[4], [0,0])
0.039275059047631206

julia> typeof(ans)
Float64

julia> Ecua[1].coeffs[4].coeffs[1]
 0.039275059047631206

julia> Ecua[1].coeffs[4].coeffs[1].coeffs[1]
0.039275059047631206

julia> typeof(Ecua[1].coeffs[4].coeffs[1])
TaylorSeries.HomogeneousPolynomial{Float64}

julia> typeof(Ecua[1].coeffs[4].coeffs[1].coeffs[1])
Float64

julia> jacobian(Ecua[1].coeffs[4])
ERROR: MethodError: `jacobian` has no method matching jacobian(::TaylorSeries.TaylorN{Float64})

julia> jacobian(Ecua[1])
ERROR: MethodError: `jacobian` has no method matching jacobian(::TaylorSeries.Taylor1{TaylorSeries.TaylorN{Float64}})

julia> jacobian(Ecua[)
ERROR: syntax: unexpected )

julia> jacobian(Ecua)
ERROR: MethodError: `jacobian` has no method matching jacobian(::Array{TaylorSeries.Taylor1{TaylorSeries.TaylorN{Float64}},1})

julia> jacobian( [Ecua[1].coeffs[4], Ecua[2].coeffs[4]] )
2x2 Array{Float64,2}:
 -3.76975   0.0    
  2.3      -5.06975






















