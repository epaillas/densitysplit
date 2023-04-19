using CellListMap
using StaticArrays
using LinearAlgebra


function _count_pairs!(i, j, d2, weights1, weights2, counts)
    if d2 > 0
        counts[i] += weights1[i] * weights2[j]
    end
    return counts
end 

function _count_pairs_gaussian!(i, j, d2, weights1, weights2, smooth_radius, counts)
    if d2 > 0
        norm = smooth_radius^3 * (2 * pi)^(3/2)
        norm = 1 
        gaussian_weight = (1/norm) * exp.(-d2 / (2 * smooth_radius^2))
        counts[i] += weights1[i] * weights2[j] * gaussian_weight
    end
    return counts
end 

function count_pairs_survey_tophat(
    positions1, positions2, weights1, weights2, smooth_radius
)
    rmax = smooth_radius
    positions1 = convert(Array{Float64}, positions1)
    positions2 = convert(Array{Float64}, positions2)
    weights1 = convert(Array{Float64}, weights1)
    weights2 = convert(Array{Float64}, weights2)
    npos1 = size(positions1)[2]
    D1D2 = zeros(Float64, npos1);
    box = Box(limits(positions1, positions2), rmax)
    cl = CellList(positions1, positions2, box)
    map_pairwise!(
        (x, y, i, j, d2, D1D2) ->
        _count_pairs!(i, j, d2, weights1, weights2, D1D2),
        D1D2, box, cl,
        parallel=true
    )
    return D1D2
end

function count_pairs_survey_gaussian(
    positions1, positions2, weights1, weights2, smooth_radius
)
    rmax = 3 * smooth_radius
    positions1 = convert(Array{Float64}, positions1)
    positions2 = convert(Array{Float64}, positions2)
    weights1 = convert(Array{Float64}, weights1)
    weights2 = convert(Array{Float64}, weights2)
    npos1 = size(positions1)[2]
    D1D2 = zeros(Float64, npos1);
    box = Box(limits(positions1, positions2), rmax)
    cl = CellList(positions1, positions2, box)
    map_pairwise!(
        (x, y, i, j, d2, D1D2) ->
        _count_pairs_gaussian!(i, j, d2, weights1, weights2, smooth_radius, D1D2),
        D1D2, box, cl,
        parallel=true
    )
    return D1D2
end


function count_pairs_box_tophat(
    positions1, positions2, weights1, weights2, box_size, smooth_radius
)
    rmax = smooth_radius
    positions1 = convert(Array{Float64}, positions1)
    positions2 = convert(Array{Float64}, positions2)
    weights1 = convert(Array{Float64}, weights1)
    weights2 = convert(Array{Float64}, weights2)
    npos1 = size(positions1)[2]
    D1D2 = zeros(Float64, npos1);
    Lbox = [box_size, box_size, box_size]
    box = Box(Lbox, rmax)
    cl = CellList(positions1, positions2, box)
    map_pairwise!(
        (x, y, i, j, d2, D1D2) ->
        _count_pairs!(i, j, d2, weights1, weights2, D1D2),
        D1D2, box, cl,
        parallel=true
    )
    return D1D2
end


function count_pairs_box_gaussian(
    positions1, positions2, weights1, weights2, box_size, smooth_radius
)
    rmax = 3 * smooth_radius
    positions1 = convert(Array{Float64}, positions1)
    positions2 = convert(Array{Float64}, positions2)
    weights1 = convert(Array{Float64}, weights1)
    weights2 = convert(Array{Float64}, weights2)
    npos1 = size(positions1)[2]
    D1D2 = zeros(Float64, npos1);
    Lbox = [box_size, box_size, box_size]
    box = Box(Lbox, rmax)
    cl = CellList(positions1, positions2, box)
    map_pairwise!(
        (x, y, i, j, d2, D1D2) ->
        _count_pairs_gaussian!(i, j, d2, weights1, weights2, smooth_radius, D1D2),
        D1D2, box, cl,
        parallel=true
    )
    return D1D2
end
