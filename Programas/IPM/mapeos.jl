module mapeos

export Estandar, EstandarI, Jung, Henon, HenonI
#-------------------------------------------------------------------------------
"""
Estandar(x::Any, y::Any k::Real)

Función para el mapeo estándar con módulo 2pi.
"""
function Estandar(x::Any,y::Any,k::Real)
    x_n = mod(x+k*sin(y),2*pi)
    y_n = mod(x+y+k*sin(y),2*pi)
    return [x_n,y_n]
end

#-------------------------------------------------------------------------------

"""
EstandarI(x::Any, y::Any k::Real)

Función para el mapeo estándar inverso con módulo 2pi.
"""
function EstandarI(x::Any, y::Any ,k::Real)
    x_n = mod(x-y+k*sin(x),2*pi)
    y_n = mod(y-k*sin(x),2*pi)
    return [x_n,y_n]
end

#-------------------------------------------------------------------------------
"""
Jung(x::Any, y::Any, a::Real )

Función para el mapeo exponencial.
"""
function Jung(x::Any , y::Any, a::Real)
    x_n = x+y
    y_n = y+a*x_n*(x_n-1.)*exp(-x_n)
    return [x_n,y_n]
end

#-------------------------------------------------------------------------------
"""
Henon(x::Any,y::Any,a::Real,b::Real)

EL mapeo de Hénon con dos parámetros.
"""
function Henon(x::Any, y::Any, a::Real, b::Real)
    x_n = a-b*y-x^2
    y_n = x
    return [x_n,y_n]
end

#-------------------------------------------------------------------------------

"""
HenonI(x::Any, y::Any,a::Real,b::Real)
Función que implementa el mapeo de Hénon inverso
"""
function HenonI(x::Any, y::Any,a::Real, b::Real)
    x_n = y
    y_n = (x+y^2-a)/(-b)
    return [x_n,y_n]
end
#-------------------------------------------------------------------------------


end
