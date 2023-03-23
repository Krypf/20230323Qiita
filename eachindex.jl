#%%
using Random
using CSV
using DataFrames
#%%
function bittoM11(a)
    2 * a - 1
end


function ising_hamiltonian1D(A; args=(-1,)) 
    J = args[1]
    _ans = J * A[end] * A[1]
    for j in eachindex(A)[1:end-1]
        _ans += J * A[j] * A[j + 1]
        # @show (j, _ans)
    end
    return _ans
end
# @time ising_hamiltonian1D(A)

#%% metropolis method
function narashi(x)# ならし
    x*(x>=0) + 0(x<0)# 非負で恒等写像, 負でゼロ写像．
end

function metropolis(s,H,β;args=(-1,),num=1, Hamiltonian=ising_hamiltonian1D)
    counts = 0
    S = deepcopy(s)
    r = 1
    p = 0
    _length = length(s)
    while (r > p)
        S = deepcopy(s)# reset needed
        j = rand(1:_length,num);
        S[j] = -s[j]
        H2 = Hamiltonian(S;args=args)# J, K are needed
        ΔE = H2 - H
        p = exp(-β * narashi(ΔE))
        r = rand()
        counts +=1

        print("\t")
        @show (j,ΔE,r,p,counts)
        if counts > 10
            println((ΔE,r,p,counts))
            break
        end
    end
    if r≤p
        #println("change")
        print("\t")
        @show (r,p)
        s = S
    end
    return s
end


#%%
M = 5000
L = 50
frame = Int(M/L)
function main(M,s,T;args=(-1,0),num=1, Hamiltonian=ising_hamiltonian1D)
    β = 1 / (T)# inverse temp. #k_B *
    for n = 1:M
        H = Hamiltonian(s;args=args); #println(H)
        @show s, H
        if (1.0 ∉ s||- 1.0 ∉ s)# all is the same
            print("all ")
            # continue
            @show n
            break
        else
            s = metropolis(s,H,β;args=args,num=num,Hamiltonian=Hamiltonian)
        end
        if n%frame==1
            q = Int((n-1)/frame)+1
            print(q, " ")
            println(s)
            # df = DataFrame(s, :auto)
            # CSV.write(string(q) * ".csv", df, writeheader=false)
        end
    end
    println(Hamiltonian(s;args=args))
    return nothing # s
end
#%%
N = 10 # 1000# 2^6
s = bittoM11.(bitrand(N))
# https://docs.julialang.org/en/v1/stdlib/Random/

#println("Enter T\n")
T = 1 #parse(Float64, readline())
J = -1; K = 0;
args = (J,K)
@time s = main(M,s,T; args=args)
@show (N, M, T, J, K)
#%%
