import click

from .experiments import Experiments


@click.command()
@click.option('-n', type=str, help='Values for # of facets')
@click.option('-d', type=str, help='Values for dimensions.')
@click.option('--num-bits', type=int, help='LSH num bits. default = d', default=0)
@click.option('-o', '--output', type=str, default=None,
              help='Location to store the results')
def experiment(n, d, num_bits, output):
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
    exp.run(n, d, num_bits, output)


if __name__ == '__main__':
    experiment()
