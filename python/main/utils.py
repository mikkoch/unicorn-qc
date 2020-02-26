import os
import hail as hl
import numpy as np
import pandas as pd

import matplotlib
import matplotlib.pyplot as plt

from typing import Union
from google.cloud import storage
from hail.expr.expressions.typed_expressions import ArrayExpression


def expand_pd_array_col(
        df: pd.DataFrame,
        array_col: str,
        num_out_cols: int = 0,
        out_cols_prefix=None,
        out_1based_indexing: bool = True
) -> pd.DataFrame:
    """
    Expands a Dataframe column containing an array into multiple columns.
    :param DataFrame df: input dataframe
    :param str array_col: Column containing the array
    :param int num_out_cols: Number of output columns. If set, only the `n_out_cols` first elements of the array column are output.
                             If <1, the number of output columns is equal to the length of the shortest array in `array_col`
    :param out_cols_prefix: Prefix for the output columns (uses `array_col` as the prefix unless set)
    :param bool out_1based_indexing: If set, the output column names indexes start at 1. Otherwise they start at 0.
    :return: dataframe with expanded columns
    :rtype: DataFrame
    """

    if out_cols_prefix is None:
        out_cols_prefix = array_col

    if num_out_cols < 1:
        num_out_cols = min([len(x) for x in df[array_col].values.tolist()])

    cols = ['{}{}'.format(out_cols_prefix, i + out_1based_indexing) for i in range(num_out_cols)]
    df[cols] = pd.DataFrame(df[array_col].values.tolist())[list(range(num_out_cols))]

    return df




def to_pandas_and_explode_array(x: hl.Table,
                                column_key_name: object = 's',
                                column_value_name: object = 'scores', new_column_value_name: object = 'PC',
                                drop_key: object = False, as_array: object = False,
                                key_by_s: object = True) -> object:
    """Convert hail Table to pandas DataFrame

        Parameters
        ----------
        x : hail Table
            Each row corresponding to one individual. The key for x is supposed to be "s" (sample id).
            x should also consist a column "score" in ArrayExpression type
            that stores principal components for each individual

        column_key_name: str (default = 's')
            Name for column key

        column_value_name: str (default = 'scores')
            Column name that corresponds to principal components

        drop_key: bool (default = False)
            if it is True, the key column will be dropped from return value

        as_array: bool (default = False)
            if it is True, return value will be a numpy array

        bey_by_s: bool
            if it is True, the index of pd.DataFrame will be `column_key`

        Returns
        -----------

        rev: a numpy array or pd.DataFrame
            Each column of the rev corresponds to sample Id or PC scores
        """
    x = x.to_pandas()
    x_key = x[[column_key_name]].values
    x_value = np.array(x[column_value_name].apply(pd.Series).values.tolist())
    n_features = x_value.shape[1]

    if as_array and drop_key:
        rev = x_value

    elif as_array and not drop_key:
        rev = np.concatenate([x_key, x_value], axis=1)

    elif not as_array and drop_key:
        rev = pd.DataFrame(x_value)
        rev.rename(axis=1, inplace=True, mapper={i:"%s%d"%(new_column_value_name, i) for i in range(n_features)})

    else:
        rev = np.concatenate([x_key, x_value], axis=1)
        rev = pd.DataFrame(rev)
        rev.rename(axis=1, inplace=True, mapper={0:column_key_name})
        rev.rename(axis=1, inplace=True, mapper={i:"%s%d"%(new_column_value_name, i) for i in range(1, n_features+1)})
        if key_by_s:
            rev.set_index(column_key_name, inplace=True)

    return rev


def discrete_cmap(x, base_cmap='Dark2'):
    """Define a discrete color map from a list or a defined N

    Parameters:
    ----------

    x: np.Series or a list or a set
        a Series or a list contains color labels
    base_cmap: name of cmap name

    Returns:
    ----------

    Return a dictionary<label, color>
    """

    if isinstance(x, set):
        x = list(x)
    elif isinstance(x, pd.Series):
        x = list(set(x.tolist()))
    elif isinstance(x, list):
        x = list(set(x))
    else:
        raise ValueError

    n = len(x)
    base = plt.cm.get_cmap(base_cmap)
    color_list = base(np.linspace(0, 1, n))
    cmap_name = base.name + str(n)
    color_map = matplotlib.colors.LinearSegmentedColormap.from_list('mcm', color_list, n)
    color_map = {x[i]:color_map(i) for i in range(n)}
    return(color_map)



def get_bucket(path: str)->storage.bucket.Bucket:
    """Return bucket object

    Parameters:
    ----------
    path: a google storage path

    Returns: storage.bucket.Bucket
    """
    if path.startswith('gs://'):
        path = path.replace('gs://', '')
    bucket_name = path.split('/').pop(0)
    try:
        client = storage.Client()
        bucket = client.get_bucket(bucket_name)
        return bucket
    except:
        raise ValueError('The bucket %s is not found' % bucket_name)


def get_blob(path: str)->storage.blob.Blob:
    """Return blob object

    Parameters:
    ----------
    path: a google storage path

    Returns: storage.blob.Blob
    """
    if path.startswith('gs://'):
        path = path.replace('gs://', '')
    path_fields = path.split('/')
    bucket_name = path_fields.pop(0)
    blob_name = '/'.join(path_fields)
    try:
        client = storage.Client()
        bucket = client.get_bucket(bucket_name)
        blob = bucket.blob(blob_name)
        return blob
    except:
        raise ValueError('The bucket %s is not found' % bucket_name)



def is_dir(blob:Union[str, storage.blob.Blob])->bool:
    """Return whether the blob links to a directory or a file
        The basic idea is that directory ends with '/'

    Parameters:
    ----------
    blob: name of the blob or blob in storage.blob.Blob type

    Returns:
    ----------
    whether it is a blob
    """
    if isinstance(blob, storage.blob.Blob):
        blob_name = blob.name
    elif isinstance(blob, str):
        blob_name = blob
    else:
        raise ValueError('The input type %s is not allowed' % type(blob))
    return blob_name.endswith('/')


def is_parent(candidate_parent_path: str, child_path: str, recursive: bool=False)->bool:
    """ Return whether is parent

    Parameters:
    ----------
    candidate_parent_path: str
    child_path: str
    Returns:
    ----------
    """
    if candidate_parent_path.startswith('gs://'):
        candidate_parent_path = candidate_parent_path.replace('gs://', '')
    if candidate_parent_path.endswith('/'):
        candidate_parent_path = candidate_parent_path[:-1]
    if child_path.startswith('gs://'):
        child_path = child_path.replace('gs://', '')
    if child_path.endswith('/'):
        child_path = child_path[:-1]
    if recursive:
        return candidate_parent_path in child_path
    else:
        return os.path.dirname(child_path) == candidate_parent_path



def list_dirs(path: str, local_run: bool=False, recursive:bool=False, add_gs_prefix:bool=True,
               credential_path='/Users/danfengc/Documents/UNICORN-671b21bc7295.json')->list:
    """List all files (not directories under a google bucket)

        Parameters:
        ----------
        google_bucket: str
            bucket name to search

        Returns:
        ----------
        A list of files that belongs to the google bucket
        """
    if local_run:
        os.environ['GOOGLE_APPLICATION_CREDENTIALS'] = credential_path

    rev = []
    bucket = get_bucket(path)
    blobs = bucket.list_blobs()

    for blob in blobs:
        if not is_dir(blob):
            continue
        if add_gs_prefix:
            filename = os.path.join('gs://', bucket.name, blob.name[:-1])
        else:
            filename = os.path.join(bucket.name, blob.name[:-1])

        if is_parent(candidate_parent_path=path, child_path=filename, recursive=recursive):
            rev.append(filename)

    return rev


def list_files(path: str, local_run: bool=False, recursive:bool=False, add_gs_prefix:bool=True,
               credential_path='/Users/danfengc/Documents/UNICORN-671b21bc7295.json')->list:
    """List all files (not directories under a google bucket)

    Parameters:
    ----------
    google_bucket: str
        bucket name to search

    Returns:
    ----------
    A list of files that belongs to the google bucket
    """
    if local_run:
        os.environ['GOOGLE_APPLICATION_CREDENTIALS'] = credential_path

    rev = []
    bucket = get_bucket(path)
    blobs = bucket.list_blobs()

    for blob in blobs:
        if is_dir(blob):
            continue
        if add_gs_prefix:
            filename = os.path.join('gs://', bucket.name, blob.name)
        else:
            filename = os.path.join(bucket.name, blob.name)
        if is_parent(candidate_parent_path=path, child_path=filename, recursive=recursive):
            rev.append(filename)

    return rev


def file_exists(file_path: str, local_run: bool=False,
                credential_path='/Users/danfengc/Documents/UNICORN-671b21bc7295.json')->bool:
    """Function for testing whether a file exists in google storage

    Parameters:
    ----------
    file_path: str
        file path of google storage

    Returns:
    ----------
    boolean value indicating whether file exists or not
    """
    if local_run:
        os.environ['GOOGLE_APPLICATION_CREDENTIALS'] = credential_path

    blob = get_blob(file_path)
    return blob.exists()
