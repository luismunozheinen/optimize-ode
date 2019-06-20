module setup_ode_def

using DifferentialEquations,EqSimRSE_Geometry,ParameterizedFunctions
export frs!

#Functions---------------------------------------------------------------------
#1. Interactions/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
function LS_connect(d,l)::Array{Int64,2}
    ls=zeros(d*l-d,2)
    i=1
    for row=1:d
        for col=1:(l-1)
        ls[i,1]=col+(row-1)*l
        ls[i,2]=col+1+(row-1)*l
        i=i+1
        end
    end
    return ls
end

function LS_connect_closed(d,l)::Array{Int64,2}
    ls=zeros(d*l,2)
    i=1
    for row=1:d
        for col=1:(l-1)
        ls[i,1]=col+(row-1)*l
        ls[i,2]=col+1+(row-1)*l
        i=i+1
        end
        if l!=1
        ls[i,1]=1+(row-1)*l
        ls[i,2]=l+(row-1)*l
        i=i+1
        end
    end
    return ls
end

function CS_connect(d,l)::Array{Int64,2}
    cs=zeros(d*l-l,2)
    i=1
    for row=1:(d-1)
        for col=1:l
        cs[i,1]=col+(row-1)*l
        cs[i,2]=col+row*l
        i=i+1
        end
    end
    return cs
end
function build_links!(connect)
    #Setup ODE String
    #lc=LS_connect(d,l)
    lc=LS_connect_closed(d,l)
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
#2. Equations/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
#Setup individual equations for du, dv and dθ
function write_du_rs(i)
    "du$i=v$i-ν\n"
end
function write_dv_rs(i)
        string("dv$i=-1/m*(kp[$i]*u$i+σ[$i]*A[$i]*asinh(v$i/(2*ν)*exp(θ$i/A[$i]))-kc[$i]*(",connect[i],"))\n")
end
function write_dθ_rs(i)
    "dθ$i=B[$i]*ν/L*(exp((τ0-θ$i)/B[$i])-v$i/ν)\n"
end
#Assemble System----------------------------------------------------------------
#Load Connectivities
connect=Array{String}(undef,d*l)
connect[:].=""
build_links!(connect)
#Assemble system in string "odeSystem"
odeSystem=""
#Load each line of the system
for i=1:d*l
    global odeSystem
    odeSystem=string(odeSystem,write_du_rs(i),write_dv_rs(i),write_dθ_rs(i))
end
writeFunction=string("frs! = @ode_def StickSlip_RSE_dim begin\n",odeSystem,"end m kp kc ν σ τ0 L A B ")
#Convert string into expression for ode system
line=Meta.parse(writeFunction)
#Evaluate frs!
eval(line)
end
