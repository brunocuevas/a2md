from fragbuilder import Peptide, set_seed
import openbabel
import click
import sys
import time

@click.group()
def cli():
    pass

@click.command()
@click.option('--output', default=None, help='file to save the info')
@click.option('--number', default=1, help='number of conformers')
@click.option('--seed', default=42, help='seed')
@click.option('--sample', default=1, help='residue in which to perform sampling')
@click.argument('sequence')
def generate_peptide(sequence, sample, output, seed, number):

    start = time.time()
    set_seed(seed)

    for i in range(number):
        pept = Peptide(sequence.upper(), nterm="neutral", cterm="neutral")
        pept.sample_bb_angles(sample)
        if sequence[sample] != "P":
            pept.sample_chi_angles(sample)
        if output is None:
            name = sequence + '_{:06d}.pdb'.format(i)
        else:
            name = output + '_{:06d}.pdb'.format(i)

        pept.write_pdb(str(name))
        bb = pept.get_bb_angles(2)
        print(bb)

    # print("TE : {:12.4f}".format(time.time() - start))
    return  True

@click.command()
@click.option('--output', default=None, help='file to save the info')
@click.option('--number', default=1, help='number of conformers')
@click.option('--seed', default=42, help='seed')
@click.option('--cterm', default='methyl', help='either methyl, charged or neutral')
@click.option('--nterm', default='methyl', help='either methyl, charged or neutral')
@click.argument('residue')
def generate_residue(residue, output, seed, number, cterm, nterm):

    set_seed(seed)

    for i in range(number):
        pept = Peptide(residue.upper(), nterm=nterm, cterm=cterm)
        pept.sample_bb_angles(1)
        if residue != "P":
            pept.sample_chi_angles(1)
        if output is None:
            name = 'c{:s}c'.format(residue) + '_{:06d}.mol2'.format(i)
        else:
            name = output + '_{:06d}.mol2'.format(i)

        pept.write_file('mol2', str(name))
        bb = pept.get_bb_angles(1)
        print(bb)

    return  True


cli.add_command(generate_peptide)
cli.add_command(generate_residue)

if __name__ == "__main__":
    cli()