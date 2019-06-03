from __future__ import print_function

import itertools
import json
from pprint import pprint
from timeit import default_timer as timer

import numpy as np

import volesti
from .utils import random_polytope


class Experiments:

    @staticmethod
    def run_single(p, n_queries=1000, num_bits=0):
        internal_point = np.zeros(p.dimension)

        tic = timer()
        print(f'asd  -- {internal_point.shape}')
        p.create_point_representation(internal_point.squeeze(), num_bits)
        print('asd2')
        toc = timer()
        constr_time = toc - tic

        queries = p.sample_boundary(internal_point, 2 * n_queries, 5)
        for idx in range(n_queries, 2 * n_queries):
            queries[idx] = queries[idx] * 1.1

        correct = 0
        avg_time = 0
        is_in_avg_time = 0
        for query in queries:
            tic = timer()
            contains_naive = p.is_in(query)
            toc = timer()
            is_in_avg_time += toc - tic

            tic = timer()
            contains_oracle = p.contains_point(query)
            toc = timer()
            avg_time += toc - tic

            if contains_naive == contains_oracle:
                correct += 1

        return constr_time, avg_time / len(queries), is_in_avg_time / len(queries), correct / len(queries)

    @staticmethod
    def on_the_fly(n_d, n_queries=1000, num_bits=0):
        n, d = n_d
        p = random_polytope(n, d)
        return Experiments.run_single(p, n_queries, num_bits)

    @staticmethod
    def offline(p, n_queries=1000, num_bits=0):
        return Experiments.run_single(p, n_queries, num_bits)

    @staticmethod
    def run(polytopes=None, num_bits=0, output=None, n_values=None, d_values=None, n_iter=5, vis=True):
        results = []
        iter_idx = 0

        if n_values is not None and d_values is not None:
            the_values = itertools.product(n_values, d_values)
            function = Experiments.on_the_fly
        else:
            the_values = polytopes
            function = Experiments.offline

        for value in the_values:
            if isinstance(value, tuple):
                n, d = value
            else:
                n, d = value.n_hyperplanes, value.dimension
            pair_avg_constr_time = 0
            pair_avg_query_time = 0
            pair_avg_is_in_time = 0
            pair_avg_success_rate = 0

            for _ in range(n_iter):
                print('{} out of {}'.format(iter_idx + 1, (len(d_values) if d_values else 1) * (len(n_values) if n_values else 1) * n_iter))
                constr_time, query_time, is_in_time, success_rate = function(value, num_bits=num_bits)
                pair_avg_constr_time += constr_time
                pair_avg_query_time += query_time
                pair_avg_is_in_time += is_in_time
                pair_avg_success_rate += success_rate
                iter_idx += 1
            pair_avg_constr_time /= n_iter
            pair_avg_query_time /= n_iter
            pair_avg_is_in_time /= n_iter
            pair_avg_success_rate /= n_iter
            results.append({
                'd': d,
                'n': n,
                'pair_constr_time': pair_avg_constr_time,
                'pair_query_time': pair_avg_query_time,
                'pair_is_in_time': pair_avg_is_in_time,
                'pair_success_rate': pair_avg_success_rate
            })
        if output:
            with open(output, 'w') as f:
                json.dump(results, f)

        if vis:
            Experiments.visualize(results)

        return results

    @staticmethod
    def visualize(results):
        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        data = [(values['d'], values['n'], values['pair_query_time']) for values in results]
        data_is_in = [(values['d'], values['n'], values['pair_is_in_time']) for values in results]
        ax.scatter([d[0] for d in data], [d[1] for d in data], [d[2] for d in data], label='Oracle queries', color='red')
        ax.scatter([d[0] for d in data_is_in], [d[1] for d in data_is_in], [d[2] for d in data_is_in], label='naive queries', color='blue')
        ax.set_xlabel('dimension')
        ax.set_ylabel('# of facets')
        ax.set_zlabel('avg query time')
        ax.legend()

        plt.show()
