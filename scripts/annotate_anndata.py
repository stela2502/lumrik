#!/usr/bin/env python3
"""Annotate an existing H5AD in place from a cluster-to-cell-type TSV.

Usage: annotate_anndata.py FILE.h5ad ANNOTATIONS.tsv
Required TSV columns: <cluster_column>, cell_type, lineage, developmental_order
Only /obs is modified. No AnnData, pandas, or matrix loading.
"""
import argparse
import csv
import sys
from pathlib import Path

FIELDS = ('cell_type', 'lineage', 'developmental_order')
STAGING = '__lumrik_annotation_staging__'

class AnnotationError(Exception):
    pass

def as_text(value):
    if isinstance(value, bytes):
        return value.decode('utf-8')
    return str(value)

def read_table(path):
    try:
        with path.open(encoding='utf-8-sig', newline='') as fh:
            reader = csv.DictReader(fh, delimiter='\t')
            cluster_column = reader.fieldnames[0] if reader.fieldnames else None
            required = (*FIELDS,)
            if not reader.fieldnames or len(reader.fieldnames) != len(set(reader.fieldnames)):
                raise AnnotationError('TSV header is missing or has duplicate column names')
            absent = [name for name in required if name not in reader.fieldnames]
            if absent:
                raise AnnotationError('TSV missing required columns: ' + ', '.join(absent))
            if cluster_column in FIELDS:
                raise AnnotationError('First TSV column must identify the H5AD cluster column, not an annotation field')
            result = {}
            for line, row in enumerate(reader, 2):
                if None in row or any(v is None for v in row.values()):
                    raise AnnotationError(f'TSV line {line}: incorrect number of columns')
                cid = row[cluster_column].strip()
                if not cid:
                    raise AnnotationError(f'TSV line {line}: empty cluster')
                if cid in result:
                    raise AnnotationError(f'TSV line {line}: duplicate cluster {cid!r}')
                cell_type, lineage, order = (row[k].strip() for k in FIELDS)
                if not cell_type or not lineage or not order:
                    raise AnnotationError(f'TSV line {line}, cluster {cid}: empty annotation field')
                try:
                    numeric_order = int(order)
                except ValueError:
                    raise AnnotationError(f'TSV line {line}, cluster {cid}: developmental_order must be an integer') from None
                if numeric_order < 0:
                    raise AnnotationError(f'TSV line {line}, cluster {cid}: developmental_order must be >= 0')
                result[cid] = (cell_type, lineage, numeric_order)
            if not result:
                raise AnnotationError('Annotation TSV contains no clusters')
            return cluster_column, result
    except OSError as exc:
        raise AnnotationError(f'Cannot read TSV {path}: {exc}') from None

def cluster_ids(obs, cluster_column):
    if cluster_column not in obs:
        raise AnnotationError(f'H5AD /obs/{cluster_column} is missing; cluster annotation requires this column')
    ds = obs[cluster_column]
    if hasattr(ds, 'keys'):
        if 'codes' not in ds or 'categories' not in ds:
            raise AnnotationError(f'Unsupported /obs/{cluster_column} categorical encoding (missing codes/categories)')
        categories = [as_text(x) for x in ds['categories'].asstr()[:]]
        codes = ds['codes'][:]
        if any(int(c) < 0 or int(c) >= len(categories) for c in codes):
            raise AnnotationError(f'/obs/{cluster_column} contains missing or invalid categorical codes')
        return [categories[int(c)] for c in codes]
    if getattr(ds.dtype, 'kind', None) in ('S', 'O', 'U'):
        return [as_text(x) for x in ds.asstr()[:]]
    return [str(x) for x in ds[:]]

def annotate(h5path, tsvpath):
    try:
        import h5py
        import numpy as np
    except (ImportError, ValueError, RuntimeError, AttributeError) as exc:
        raise AnnotationError(f'Cannot import h5py/NumPy: {type(exc).__name__}: {exc}. Use a compatible Python environment.') from None
    if not h5path.is_file():
        raise AnnotationError(f'H5AD not found: {h5path}')
    if h5path.resolve() == tsvpath.resolve():
        raise AnnotationError('H5AD and annotation TSV must be different files')
    cluster_column, mapping = read_table(tsvpath)
    try:
        with h5py.File(h5path, 'r+') as h5:
            if 'obs' not in h5:
                raise AnnotationError('H5AD has no /obs group')
            obs = h5['obs']
            ids = cluster_ids(obs, cluster_column)
            actual = set(ids)
            missing = sorted(actual - mapping.keys())
            unknown = sorted(mapping.keys() - actual)
            if missing or unknown:
                raise AnnotationError(f'clusters do not match TSV. Missing: {missing}; unknown: {unknown}')
            if STAGING in obs:
                raise AnnotationError(f'/obs/{STAGING} exists from a previous interrupted run; inspect before retrying')
            raw_order = obs.attrs.get('column-order')
            if raw_order is None:
                raise AnnotationError('/obs lacks AnnData column-order metadata; refusing to modify unknown layout')
            order = [as_text(x) for x in raw_order]
            if any((name in obs) != (name in order) for name in FIELDS):
                raise AnnotationError('Existing annotation datasets disagree with /obs column-order metadata')
            if len(ids) != len(obs[as_text(obs.attrs['_index'])]):
                raise AnnotationError(f'/obs/{cluster_column} length does not match obs index length')
            # Stage only small annotation datasets. Never touch /X or copy the H5AD.
            staging = obs.create_group(STAGING)
            old_names = [name for name in FIELDS if name in obs]
            published = []
            try:
                utf8 = h5py.string_dtype(encoding='utf-8')
                for pos, name in enumerate(FIELDS):
                    values = [mapping[cid][pos] for cid in ids]
                    if name == 'developmental_order':
                        ds = staging.create_dataset(name, data=np.asarray(values, dtype=np.int64))
                        ds.attrs['encoding-type'] = 'array'
                        ds.attrs['encoding-version'] = '0.2.0'
                    else:
                        # A factor is a standard AnnData categorical: codes + categories.
                        # Cell types use developmental_order for their level order.
                        if name == 'cell_type':
                            first_order = {}
                            for cell_type, _lineage, dev_order in mapping.values():
                                first_order[cell_type] = min(first_order.get(cell_type, dev_order), dev_order)
                            levels = sorted(set(values), key=lambda v: (first_order[v], v))
                            ordered = True
                        else:
                            levels = sorted(set(values))
                            ordered = False
                        codes_by_level = {level: i for i, level in enumerate(levels)}
                        grp = staging.create_group(name)
                        grp.attrs['encoding-type'] = 'categorical'
                        grp.attrs['encoding-version'] = '0.2.0'
                        grp.attrs['ordered'] = np.bool_(ordered)
                        grp.create_dataset('categories', data=np.asarray(levels, dtype=object), dtype=utf8)
                        grp['categories'].attrs['encoding-type'] = 'string-array'
                        grp['categories'].attrs['encoding-version'] = '0.2.0'
                        grp.create_dataset('codes', data=np.asarray([codes_by_level[v] for v in values], dtype=np.int32))
                        grp['codes'].attrs['encoding-type'] = 'array'
                        grp['codes'].attrs['encoding-version'] = '0.2.0'
                h5.flush()
                # Move existing columns into staging as backups before publishing.
                for name in old_names:
                    obs.move(name, f'{STAGING}/__old_{name}')
                for name in FIELDS:
                    obs.move(f'{STAGING}/{name}', name)
                    published.append(name)
                obs.attrs['column-order'] = np.asarray(
                    order + [name for name in FIELDS if name not in order], dtype=utf8
                )
                h5.flush()
            except Exception:
                for name in published:
                    if name in obs:
                        del obs[name]
                for name in old_names:
                    backup = f'{STAGING}/__old_{name}'
                    if backup in obs:
                        obs.move(backup, name)
                obs.attrs['column-order'] = np.asarray(order, dtype=utf8)
                if STAGING in obs:
                    del obs[STAGING]
                h5.flush()
                raise
            # All new columns are committed. Remove small obsolete annotations.
            del obs[STAGING]
            h5.flush()
            print(f'Annotated {h5path}: {len(ids):,} cells, {len(actual)} clusters; updated /obs only.')
    except OSError as exc:
        raise AnnotationError(f'Cannot update H5AD {h5path}: {exc}') from None

def main():
    parser = argparse.ArgumentParser(description='Annotate an H5AD in place using a cluster-to-cell-type TSV. Never copies the H5AD.')
    parser.add_argument('h5ad', type=Path, help='Existing H5AD (modified in place)')
    parser.add_argument('annotations', type=Path, help='TSV columns: <cluster_column>, cell_type, lineage, developmental_order')
    args = parser.parse_args()
    try:
        annotate(args.h5ad, args.annotations)
    except AnnotationError as exc:
        print(f'ERROR: {exc}', file=sys.stderr)
        return 1
    except Exception as exc:
        print(f'ERROR: unexpected {type(exc).__name__}: {exc}', file=sys.stderr)
        return 1
    return 0

if __name__ == '__main__':
    sys.exit(main())
