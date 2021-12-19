from os import path
import numpy as np
from julia.api import Julia


def generate_seeds(
    nseeds, box_size=None, selection_function='uniform',
    sampling_data=None
):
    if selection_function == 'randoms':
        # sample from randoms file
        idx = np.random.choice(
            len(sampling_data), size=nseeds, replace=False
        )
        seeds = sampling_data[idx]
    else:
        # sample from a uniform distribution
        x = np.random.uniform(0, box_size, nseeds)   
        y = np.random.uniform(0, box_size, nseeds)   
        z = np.random.uniform(0, box_size, nseeds) 
        seeds = np.c_[x, y, z]
    return seeds


def get_density_pdf(
    smoothing_radius, data_positions1, data_weights1,
    data_positions2, data_weights2, selection_function='uniform',
    randoms_positions2=None, randoms_weights2=None, box_size=None
):
    """
    Split the random seeds according to the local
    galaxy density.
    """

    # import Julia modules
    jl = Julia(compiled_modules=False)
    from julia import Main

    module_path = path.join(path.dirname(__file__),
        'fastmodules', 'count_pairs.jl')

    jl.eval(f'include("{module_path}")')

    if selection_function == 'randoms':
        Main.positions1 = data_positions1
        Main.weights1 = data_weights1
        Main.positions2 = data_positions2
        Main.weights2 = data_weights2
        Main.rmax = smoothing_radius

        D1D2 = jl.eval(
            "count_pairs("
            "positions1, positions2, weights1,"
            "weights2", "rmax)"
        )

        Main.positions2 = randoms_positions2 
        Main.weights2 = randoms_weights2 

        D1R2 = jl.eval(
            "count_pairs("
            "positions1, positions2, weights1,"
            "weights2", "rmax)"
        )

    else:
        Main.positions1 = data_positions1
        Main.weights1 = data_weights1 
        Main.positions2 = data_positions2
        Main.weights2 = data_weights2
        Main.box_size = box_size
        Main.rmax = smoothing_radius

        D1D2 = jl.eval(
            "count_pairs("
            "positions1, positions2, weights1,"
            "weights2", "box_size", "rmax)"
        )

        bin_volume = 4/3 * np.pi * smoothing_radius ** 3
        mean_density = len(data_positions2) / (box_size ** 3)

        D1R2 = bin_volume * mean_density * len(data_positions1)

    delta = D1D2 / D1R2 - 1

    return delta


def generate_quantiles(seeds, density_pdf, nquantiles):
    nseeds = len(seeds)
    idx = np.argsort(density_pdf)
    sorted_seeds = seeds[idx]

    quantiles = {}
    for i in range(1, nquantiles + 1):
        quantiles[f'DS{i}'] = sorted_seeds[
            int((i - 1) * nseeds / nquantiles):int(i * nseeds / nquantiles)
        ]

    return quantiles

