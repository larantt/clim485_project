import pandas as pd

GAS_ABS = ['freq_GHz', 'dry_air_tau', 'wv_abs']
LIQ_CLD = ['freq_GHz', 'mass_abs']

def read_text_file(filepath,col_names=['freq_GHz', 'dry_air_tau', 'wv_abs']):
    """
    Reads text file for absorbption coefficient or optical depths
    
    Parameters
    ----------
    filepath : str
        path to the file to read in
    col_names : list(str)
        names of the columns in the file

    Returns
    -------
    df : pd.DataFrame
        pandas dataframe containing text file data
    """
    df = pd.read_csv(
        filepath,
        delim_whitespace=True,
        comment='#',
        header=None,
        names=col_names
    )
    return df
