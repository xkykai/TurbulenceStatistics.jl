using MAT
using CairoMakie
using Statistics
using StatsBase
using JLD2

data = matread("data-1/data/Data19.mat")
keys(data)

times = vec(data["t"])
zs = vec(data["z"])
xs = vec(data["x"])
y0s = vec(data["y0"])
us = data["u"]

#%%
# function find_min(a...)
#     return minimum(minimum.([a...]))
# end

# function find_max(a...)
#     return maximum(maximum.([a...]))
# end

# frame = 5

# clim = (minimum(us[:, :, :, frame]), maximum(us[:, :, :, frame]))

# #%%
# fig = Figure(resolution=(1000, 1000))
# ax1 = Axis(fig[1, 1], xlabel="x", ylabel="z", title="y = $(y0s[1])")
# hm_1 = heatmap!(ax1, xs, zs, us[:, 1, :, frame], colorrange=clim)
# ax2 = Axis(fig[2, 1], xlabel="x", ylabel="z", title="y = $(y0s[2])")
# hm_2 = heatmap!(ax2, xs, zs, us[:, 2, :, frame], colorrange=clim)
# Colorbar(fig[1:2, 2], hm_1)
# display(fig)
#%%
ubars = mean(us, dims=(1, 3, 4))
u′s = us .- ubars
u′²bars = mean(u′s .^ 2, dims=(1, 3, 4))

function compute_corr(u′, u′²bars, ix, iz)
    Nx = size(u′)[1]
    Nz = size(u′)[3]

    i_displaced = mod1.(collect(1:Nx) .+ ix, Nx)
    k_displaced = mod1.(collect(1:Nz) .+ iz, Nz)

    return mean(u′ .* @view(u′[i_displaced, :, k_displaced, :]), dims=(1, 3, 4)) ./ u′²bars
end

function compute_corr(u′, u′²bars)
    Nx = size(u′)[1]
    Ny = size(u′)[2]
    Nz = size(u′)[3]
    ixs = 1:floor(Int(ceil(Nx / 2)))
    izs = 1:floor(Int(ceil(Nz / 2)))
    Rs = zeros(length(ixs), Ny, length(izs))

    @inbounds for ix in ixs
        @info ix
        @inbounds Threads.@threads for iz in izs
            Rs[ix, :, iz] .= compute_corr(u′, u′²bars, ix, iz)[1, :, 1, 1]
        end
    end
    return Rs, ixs, izs
end

R1 = compute_corr(u′s, u′²bars, 1, 0)
R2 = compute_corr(u′s, u′²bars, 1535, 0)

Rs, ixs, izs = compute_corr(u′s, u′²bars)

jldsave("R_data.jld2"; Rs, ixs, izs)

# file = jldopen("R_data.jld2", "r")
# close(file)
#%%