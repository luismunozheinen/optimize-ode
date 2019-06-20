module setup_A_matrix
export friction,friction!,create_matrix
#Model Size
l = 20 #Length
d = 10 #Depth
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
        ls[i,2] = col + (row-1)*l
        i = i + 1
        end
    end
    return ls
end
"""
Function creates an array including the index of blocks which are linked (horizontal) by elastic springs
Updated Boundary Conditions
"""
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
Function assembles the linear components (elastic interaction) inbetween each block into a matrix A
"""
function create_matrix(p)
    #Retrieve Block Interactions
    lc = LS_connect_closed(d,l)
    cc = CS_connect(d,l)
    #Setup Matrix
    dof = 3*l*d                 #3 degrees of freedom per block i.e. u,v,θ
    A = zeros(dof,dof)          #Initialization
    #du=v
    for i = 1:l*d
        A[3*i-2,3*i-1] = 1.0
    end
    #boundary
    top = collect(1:l) #top row
    bottom = collect((1+l*(d-1)):d*l) #last row
    left = collect(1:l:d*l) #...
    right = collect(l:l:d*l)
    #Basic Equation x(i,j) relation
    for i = 1:l*d
        A[3*i-1,3*i-2] = (-4*p[3]-p[2])/p[1]
    end
    #Change ODE on boundaries
    [A[3*i-1,3*i-2] = (-3*p[3]-p[2])/p[1] for i in vcat(top,bottom)]
    #Corners (NOT Applicable if Closed links)
    #[A[3*i-1,3*i-2] = (-2*p[3]-p[2])/p[1] for i in vcat(1,l,l*d-l,l*d)]
    #Interactions
    for i = 1:size(lc,1)
        #Horizontal Assembly
        A[3*lc[i,1]-1,3*lc[i,2]-2] = A[3*lc[i,1]-1,3*lc[i,2]-2] + p[3]/p[1]
        A[3*lc[i,2]-1,3*lc[i,1]-2] = A[3*lc[i,2]-1,3*lc[i,1]-2] + p[3]/p[1]
    end
    for i = 1:size(cc,1)
        #Vertical Assembly
        A[3*cc[i,1]-1,3*cc[i,2]-2] = A[3*cc[i,1]-1,3*cc[i,2]-2] + p[3]/p[1]
        A[3*cc[i,2]-1,3*cc[i,1]-2] = A[3*cc[i,2]-1,3*cc[i,1]-2] + p[3]/p[1]
    end
      return A
end

end
