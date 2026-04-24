import numpy as np

def read_stat_binary_file(filename):
    """
    Read a binary file with the statistics of the input signals.

    These binary files contain the following structure:
    data_format = [('magic_word', '<Q'),
                   ('mean_chA', '<f'),
                   ('mean_chB', '<f'),
                   ('mean_chC', '<f'),
                   ('mean_chD', '<f'),
                   ('std_chA', '<f'),
                   ('std_chB', '<f'),
                   ('std_chC', '<f'),
                   ('std_chD', '<f'),
                   ('min_chA', '<f'),
                   ('min_chB', '<f'),
                   ('min_chC', '<f'),
                   ('min_chD', '<f'),
                   ('max_chA', '<f'),
                   ('max_chB', '<f'),
                   ('max_chC', '<f'),
                   ('max_chD', '<f'),
                   ('time_stamp', '<Q')] # The time is in microseconds

    Parameters
    ----------
    filename: str
        Name of the binary file to read.

    Returns
    -------
    dict()
        Dictionary containing the fields present in the binary file.
    """
    data_format = [('magic_word', '<Q'),
                   ('mean_chA', '<f'),
                   ('mean_chB', '<f'),
                   ('mean_chC', '<f'),
                   ('mean_chD', '<f'),
                   ('std_chA', '<f'),
                   ('std_chB', '<f'),
                   ('std_chC', '<f'),
                   ('std_chD', '<f'),
                   ('min_chA', '<f'),
                   ('min_chB', '<f'),
                   ('min_chC', '<f'),
                   ('min_chD', '<f'),
                   ('max_chA', '<f'),
                   ('max_chB', '<f'),
                   ('max_chC', '<f'),
                   ('max_chD', '<f'),
                   ('time_stamp', '<Q')]

    output_data = dict()
    with open(filename, 'rb') as f:
        chunk = np.fromfile(f, dtype=data_format)
        # These are the fields that will be stored in the output dictionary:
        for field in ["std_chA", "std_chB", "std_chC", "std_chD", "time_stamp"]:
            output_data[field] = chunk[field]
    return output_data