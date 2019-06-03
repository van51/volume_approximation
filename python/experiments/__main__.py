import glob
import os
from collections import defaultdict
from functools import partial
from multiprocessing import Pool

import click
import itertools
import random
import numpy as np
import tqdm

from .experiments import Experiments
from .utils import polytope_dump as dump
from .utils import polytope_load as load
from .utils import random_polytope


@click.group()
def experiment():
    pass


@experiment.command()
@click.option('-n', type=str, help='Values for # of facets')
@click.option('-d', type=str, help='Values for dimensions.')
@click.option('--num-bits', type=int, help='LSH num bits. default = d', default=0)
@click.option('-o', '--output', type=str, default=None,
              help='Location to store the results')
@click.option('--seed', type=int, default=42)
def online(n, d, num_bits, output, seed):
    random.seed(seed)
    np.random.seed(seed)

    if output is not None:
        try:
            with open(output, 'w') as _:
                pass
        except IOError as e:
            print(e)
            return
    exp = Experiments()
    n = [int(value) for value in n.split(',')]
    d = [int(value) for value in d.split(',')]
    exp.run(num_bits=num_bits, output=output, n_values=n, d_values=d, n_iter=5)


def process(fname, num_bits, output):
    p = load(fname)
    return Experiments.run(num_bits=num_bits, output=output, polytopes=[p], n_iter=1, vis=False)[0]


@experiment.command()
@click.option('-p', '--pattern')
@click.option('--num-bits', type=int, help='LSH num bits. default = d', default=0)
@click.option('-o', '--output', type=str, default=None,
              help='Location to store the results')
@click.option('--seed', type=int, default=42)
@click.option('-n', '--n-threads', type=int, default=1)
def offline(pattern, num_bits, output, seed, n_threads):
    random.seed(seed)
    np.random.seed(seed)

    filenames = glob.glob(pattern)

    _process = partial(process, num_bits=num_bits, output=output)
    pool = Pool(n_threads)
    grouped_results = defaultdict(list)
    for result in pool.map(_process, filenames):
        grouped_results[(result['n'], result['d'])].append(result)

    final_results = {}
    for pair_n_d, results in grouped_results.items():
        values = defaultdict(list)
        for key, value in results.items():
            if key == 'n' or key == 'd':
                continue
            values[key].append(value)

        final_results[pair_n_d] = {
            k: sum(v) / len(v) for k, v in values.items()
        }
    print(final_results)


@experiment.command()
@click.option('-n', type=str, help='Values for # of facets')
@click.option('-d', type=str, help='Values for dimensions.')
@click.option('-o', '--output', type=str, default=None,
              help='Location to store the results')
@click.option('-i', '--n-iter', default=3, type=int, help='How many polytopes per n x d combination to generate')
@click.option('--seed', type=int, default=42)
def create_polytopes(n, d, output, n_iter, seed):
    random.seed(seed)
    np.random.seed(seed)

    n = [int(value) for value in n.split(',')]
    d = [int(value) for value in d.split(',')]

    for n_facets, dim, i in tqdm.tqdm(itertools.product(n ,d, range(n_iter)), total=len(n) * len(d) * n_iter):
        fname = os.path.join(output, 'random_polytope_{}_{}_{}'.format(n_facets, dim, i))
        dump(random_polytope(n_facets, dim), fname)


if __name__ == '__main__':
    experiment()
