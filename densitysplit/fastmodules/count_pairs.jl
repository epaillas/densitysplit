using CellListMap
using StaticArrays
using LinearAlgebra


function _count_pairs!(i, j, d2, weights1, weights2, counts)
    if d2 > 0
        counts[i] += weights1[i] * weights2[j]
    end
    return counts
end 

function count_pairs_survey(
    positions1, positions2, weights1, weights2, rmax
)
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

function count_pairs_box(
    positions1, positions2, weights1, weights2, box_size, rmax
)
    positions1 = convert(Array{Float64}, positions1)
    positions2 = convert(Array{Float64}, positions2)
    weights1 = convert(Array{Float64}, weights1)
    weights2 = convert(Array{Float64}, weights2)
    npos1 = size(positions1)[2]
    D1D2 = zeros(Int, npos1);
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
