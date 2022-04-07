#!/usr/bin/env python


import argparse
import numpy as np
import scanpy as sc


def calculate_label_proportion(
    df,
    proportion_col,
    current_label
):
    return (
        df[df[proportion_col] == current_label].shape[0] / df.shape[0]
    )


def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="""
            Generate per experiment covariates to use in downstream model.
            """
    )

    parser.add_argument(
        '-h5', '--h5_anndata',
        action='store',
        dest='h5',
        required=True,
        help='H5 AnnData file.'
    )

    parser.add_argument(
        '-e_id', '--experiment_id',
        action='store',
        dest='exp_id',
        default='experiment_id',
        help='Experiment id column in anndata obs slot. (default: %(default)s)'
    )

    parser.add_argument(
        '-prop_col', '--proportion_covariate_column',
        action='store',
        dest='cell_label',
        default='cluster',
        help='Column to calculate cell type proportions. (default: %(default)s)'
    )

    parser.add_argument(
        '-o', '--out_file',
        action='store',
        dest='out_file',
        default='adata.h5ad',
        help='Output anndata file. (default: %(default)s)'
    )

    options = parser.parse_args()

    # Load the AnnData file
    adata = sc.read_h5ad(filename=options.h5)

    label_dict = {}
    for label in adata.obs[options.cell_label].unique():
        print(
            'Calculating proportion of {} cells in each sample...'.format(label)
        )
        proportion = adata.obs.groupby(by=options.exp_id).apply(
            calculate_label_proportion,
            proportion_col=options.cell_label,
            current_label=label
        ).to_dict()

        prop_key = '{}_proportion__autogen'.format(label.replace(' ', '_'))
        adata.obs[prop_key] = [proportion[x] for x in adata.obs[options.exp_id]]

        label_dict[label] = prop_key

    # Update ann data with keys to use downstream
    adata.obs['proportion_covariate_keys__autogen'] = [
        label_dict[x] for x in adata.obs[options.cell_label]
    ]

    adata.write(
        options.out_file,
        compression='gzip'
    )


if __name__ == '__main__':
    main()
