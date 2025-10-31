import itertools
import multiprocessing
import numpy as np
import pandas as pd
from tqdm.auto import tqdm


def chunk_iterable(data, SIZE):
    """
    Generator function that chunk an iterable into lists of length SIZE

    Args
    - data: iterable
    - SIZE: int

    Returns: generator
    """
    it = iter(data)
    while True:
        x = list(itertools.islice(it, SIZE))
        if len(x) == 0:
            return
        else:
            yield x


def chunk_dict(data, SIZE):
    """
    Chunk a dictionary into smaller dictionaries of length SIZE.
    Source: https://stackoverflow.com/a/22878842
    """
    it = iter(data)
    for i in range(0, len(data), SIZE):
        yield {k: data[k] for k in itertools.islice(it, SIZE)}


def label_counts_to_matrix(label_counts, labels, dtype=None, fillna=0):
    """
    Args
    - label_counts: sequence of dict, or dict of dict. len=m
        Sequence of dict: [{label1: value, label2: value, ...},
                           {label1: value, label2: value, ...},
                           ...]
        Dict of dict: {id1: {label1: value, label2: value, ...},
                       id2: {label1: value, label2: value, ...},
                       ...}
        Note that not all inner dicts need to have the same keys,
        but all keys should be given as the second argument, `labels`.
    - labels: list-like. len=n
        Labels for columns of the returned matrix
    - dtype: dtype. default=None
        Data type of values in label_counts. Passed onto pandas.DataFrame constructor.
    - fillna: numeric. default=0
        Value to use when an entry (row) is missing a label.

    Returns: np.ndarray, shape=(m, n)
    - Rows = ids
    - Columns = labels
    - Values = values from label_counts
    """
    if isinstance(label_counts, dict):
        df = pd.DataFrame(label_counts, index=labels, dtype=dtype).T
    else:
        df = pd.DataFrame(label_counts, columns=labels, dtype=dtype)
    mat = df.fillna(fillna).values
    return mat


def label_counts_to_interaction(
    label_counts,
    labels,
    proportion=0.1,
    min_oligos=2,
    downweighting="none",
    dtype_in=None,
    dtype_out=np.uint32,
    **kwargs,
):
    """
    Count pairwise interactions between labels.

    Args
    - label_counts: iterable of dict, or dict of dict
        Each entry (inner dict) maps from a label to its count.
    - labels: sequence. len=n
        Labels. Corresponds to order of columns in returned matrix.
    - proportion: float. default=0.1
        Value between 0 and 1. Minimum proportion of total counts
        that a label must represent within an entry to be considered
        towards an interaction.
    - min_oligos: int. default=2
        Minimum count that a label must have within an entry to be
        considered towards an interaction.
    - downweighting: str. default='none'
        Downweighting strategy to estimate pairwise interactions. Adapted from
        SPRITE pipeline options.
        - 'none': no downweighting
        - 'n_over_two': for each cluster containing both target i and target j,
            increment M_ij and M_ji each by 2/n, where n is the number of targets
            in that cluster. M_ii is the number of clusters containing target i.
    - dtype_int: dtype. default=None
        Data type of values in label_counts. See dtype in label_counts_to_matrix().
    - dtype_out: dtype. default=None
        Data type of values in returned matrix.
        If None, defaults to np.uint32 if downweighting == 'none', else float.
    - **kwargs
        Passed onto label_counts_to_matrix()
        - fillna

    Returns: np.ndarray. shape=(n, n). dtype=dtype_out
      The value at index (i, j) gives the number of entries in `label_counts`
      in which labels i and j were both present with sufficient counts
      meeting the `min_oligos` and `proportion` criteria.
    """
    assert proportion >= 0 and proportion <= 1
    assert downweighting in ("none", "n_over_two")
    if dtype_out is None:
        dtype_out = np.uint32 if downweighting == "none" else float
    mat = label_counts_to_matrix(label_counts, labels, dtype=dtype_in, **kwargs)
    mask_proportion = (
        np.divide(
            mat,
            mat.sum(axis=1, keepdims=True),
            out=np.zeros_like(mat, dtype=float),
            where=mat != 0,
        )
        >= proportion
    )
    mat_thresh = ((mat >= min_oligos) & mask_proportion).astype(dtype_out)
    if downweighting == "none":
        return mat_thresh.T @ mat_thresh
    else:  # downweighting == 'n_over_two'
        with np.errstate(divide="ignore", invalid="ignore"):
            mat_thresh_weighted = mat_thresh / mat_thresh.sum(axis=1, keepdims=True) * 2
            np.nan_to_num(mat_thresh_weighted, copy=False)
            return mat_thresh.T @ mat_thresh_weighted


def label_counts_to_interaction_chunked(
    label_counts, labels, chunksize=5000, n_processes=1, **kwargs
):
    """
    Performant wrapper around label_counts_to_interaction().

    Args
    - label_counts: iterable of dict, or dict of dict.
    - labels: sequence. len=n
    - chunksize: int. default=5000
        Performance parameter; only affects runtime but not returned value.
        Number of entries within which to consider interactions at the
        same time. Size of chunks to break `label_counts` into. Limit
        based on the amount of available memory.
    - n_processes: int. default=1
        Number of parallel processes to use.
    - **kwargs
        Passed onto label_counts_to_interaction()
        - proportion
        - min_oligos
        - downweighting
        - dtype_in
        - dtype_out
        - fillna

    Returns: np.ndarray. shape=(n, n). dtype=np.uint32
      See label_counts_to_interaction().
    """
    chunks = chunk_dict if isinstance(label_counts, dict) else chunk_iterable
    try:
        n_chunks = int(np.ceil(len(label_counts) / chunksize))
    except TypeError:
        # label_counts is an interable without a __len__() method
        n_chunks = None
    pbar = tqdm(total=n_chunks)

    def update(*_):
        pbar.update()

    if n_processes > 1:
        results = []
        with multiprocessing.Pool(n_processes) as pool:
            for chunk in chunks(label_counts, chunksize):
                results.append(
                    pool.apply_async(
                        label_counts_to_interaction, (chunk, labels), kwargs, callback=update
                    )
                )
            pool.close()
            pool.join()
        pbar.close()
        return sum(result.get() for result in results)
    out = None
    for chunk in chunks(label_counts, chunksize):
        if out is None:
            out = label_counts_to_interaction(chunk, labels, **kwargs)
        else:
            out += label_counts_to_interaction(chunk, labels, **kwargs)
        pbar.update()
    pbar.close()
    return out
