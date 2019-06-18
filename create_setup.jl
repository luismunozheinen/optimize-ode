module create_setup
using ParameterizedFunctions
export f
#------------------------------------------------------------------------------------------------------
"""
Function creates an array including the index of blocks which are linked (horizontal) by elastic springs
"""
function LS_connect(d,l)::Array{Int64,2}
    ls = zeros(d*l-d,2)
    i = 1
    for row = 1:d
        for col = 1:(l-1)
        ls[i,1] = col + (row-1)*l
        ls[i,2] = col + (row-1)*l + 1
        i = i + 1
        end
    end
    return ls
end
#------------------------------------------------------------------------------------------------------
"""
Function creates an array including the index of blocks which are linked (vertical) by elastic springs
"""
function CS_connect(d,l)::Array{Int64,2}
    cs=zeros(d*l-l,2)
    i = 1
    for row = 1:(d-1)
        for col = 1:l
        cs[i,1] = col + (row-1)*l
        cs[i,2] = col + row*l
        i = i + 1
        end
    end
    return cs
end
#------------------------------------------------------------------------------------------------------
"""
Function creates litteral code specifying the linear connection between blocks and stores them in an array "connect"
"""
function build_links!(connect;d=10,l=10)
    lc=LS_connect(d,l) #load arrays containing the index
    cc=CS_connect(d,l)
    for i=1:size(lc,1)                                  #Bloc(i+-1,j+-1)
        connect[lc[i,1]]=string(connect[lc[i,1]]," + u",lc[i,2]," - u",lc[i,1])
        connect[lc[i,2]]=string(connect[lc[i,2]]," + u",lc[i,1]," - u",lc[i,2])
    end
    for i=1:size(cc,1)                                  #Bloc(i+-1,j+-1)
        connect[cc[i,1]]=string(connect[cc[i,1]]," + u",cc[i,2]," - u",cc[i,1])
        connect[cc[i,2]]=string(connect[cc[i,2]]," + u",cc[i,1]," - u",cc[i,2])
    end
end
#------------------------------------------------------------------------------------------------------
"""
Function creates litteral code specifying the first ode for each block (defines its position du)
"""
write_du(i) = "du$i = v$i - ν\n"
#------------------------------------------------------------------------------------------------------
"""
Function creates litteral code specifying the second ode for each block (defines its velocity dv)
"""
write_dv(i)=string("dv$i = -1/m * (kp * u$i + σ*A * asinh(v$i /(2*ν) * exp(θ$i/A))-kc*(",connect[i],"))\n")
#------------------------------------------------------------------------------------------------------
"""
Function creates litteral code specifying the third ode for each block (defines its frictional state dθ)
"""
write_dθ(i)="dθ$i = B * ν/L * (exp((τ0 - θ$i) /B) - v$i/ν)\n"

#------------------------------------------------------------------------------------------------------

# Main code to setup the differential equation system in f using ParameterizedFunctions
d = 10 #number of blocks (vertical)
l = 10 #number of blocks (horizontal)
#1. create litteral expression for block interactions "connect"
connect = Array{String}(undef,d*l)
connect[:] .= ""
build_links!(connect)
#2. Assemble all individual ode's into a system "odeSystem"
odeSystem = ""

for i = 1:d*l
    global odeSystem
    odeSystem = string(odeSystem,write_du(i),write_dv(i),write_dθ(i))
end
writeFunction = string("f = @ode_def StickSlip_ODE begin\n",odeSystem,"end m kp kc ν σ τ0 L A B ") #litteral expression of the ode function
eval(Meta.parse(writeFunction)) #convert into a function
end
