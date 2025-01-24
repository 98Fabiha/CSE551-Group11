import pandas as pd
from math import factorial
from itertools import combinations
import argparse
import os
import logging


def find_msu(dataframe, groups, aggregations, att):
    """
    Find and score each Minimal Sample Unique (MSU) within the dataframe
    for the specified groups
    :param dataframe: the complete dataframe of data to score
    :param groups: an array of arrays for each group of columns to test for uniqueness
    :param aggregations: an array of aggregation methods to use for the results
    :param att: the total number of attributes (QIDs) in the dataset
    :return:
    """
    df_copy = dataframe
    # 'nple' as we may be testing a group that's a single, a tuple, triple etc
    df_updates = []
    for nple in groups:
        nple = list(nple)
        cols = nple.copy()

        # Calculate the unique value counts (fK)
        cols.append('fK')
        value_counts = df_copy[nple].groupby(nple, sort=False).size()

        if 1 in value_counts.values:
            df_value_counts = pd.DataFrame(value_counts)
            df_value_counts = df_value_counts.reset_index()
            # Change the column names
            df_value_counts.columns = cols

            # Add values for fM, MSU and SUDA
            df_value_counts['fM'] = 0
            df_value_counts['suda'] = 0
            df_value_counts.loc[df_value_counts['fK'] == 1, ['fM', 'msu', 'suda']] = \
                [1, len(nple), factorial(att - len(nple))]

            # Collect results
            df_update = pd.merge(df_copy, df_value_counts, on=nple, how='left')
            df_updates.append(df_update)

    # Return results
    if len(df_updates) > 0:
        df_updates = pd.concat(df_updates)
    return df_updates


def suda(dataframe, max_msu=2, dis=0.1, columns=None):
    """
    Special Uniqueness Detection Algorithm (SUDA)
    :param dataframe:
    :param max_msu:
    :param dis:
    :param columns: the set of columns to apply SUDA to. Defaults to None (all columns)
    :return:
    """
    logger = logging.getLogger("suda")
    logging.basicConfig()

    # Get the set of columns
    if columns is None:
        columns = dataframe.columns

    for col in columns:
        if dataframe[col].nunique() < 600:
            dataframe[col] = dataframe[col].astype(pd.CategoricalDtype(ordered=True))

    att = len(columns)
    if att > 20:
        logger.warning("More than 20 columns presented; setting ATT to max of 20")
        att = 20

    # Construct the aggregation array
    aggregations = {'msu': 'min', 'suda': 'sum', 'fK': 'min', 'fM': 'sum'}
    for column in dataframe.columns:
        aggregations[column] = 'max'

    results = []
    for i in range(1, max_msu+1):
        groups = list(combinations(columns, i))
        result = (find_msu(dataframe, groups, aggregations, att))
        if len(result) != 0:
            results.append(result)

    if len(results) == 0:
        logger.info("No special uniques found")
        dataframe["suda"] = 0
        dataframe["msu"] = None
        dataframe['fK'] = None
        dataframe['fM'] = None
        return dataframe

    # Domain completion
    for result in results:
        if 'fM' not in result.columns:
            result['fM'] = 0
            result['suda'] = 0
    dataframe['fM'] = 0
    dataframe['suda'] = 0

    # Concatenate all results
    results.append(dataframe)
    results = pd.concat(results).groupby(level=0).agg(aggregations)

    results['dis-suda'] = 0
    dis_value = dis / results.suda.sum()
    results.loc[dataframe['suda'] > 0, 'dis-suda'] = results.suda * dis_value

    results = results.fillna(0)
    return results


def main():
    logger = logging.getLogger("suda")

    argparser = argparse.ArgumentParser(description='Special Uniques Detection Algorithm (SUDA) for Python.')
    argparser.add_argument('input_path', metavar='<input>', type=str, nargs=1, default='20_dataset.csv',
                           help='The name of the CSV data file to process')
    argparser.add_argument('output_path', metavar='<output>', type=str, nargs='?', default='results.csv',
                           help='The output file name')
    argparser.add_argument('m', metavar='<m>', type=int, nargs='?', default=2,
                           help='The largest minimum sample uniqueness (MSU) to test for.')
    argparser.add_argument('d', metavar='<d>', type=float, nargs='?', default=0.1,
                           help='The file-level disclosure intrusion score (DIS)')
    argparser.add_argument('c', metavar='<c>', type=str, nargs='*', default=None, action='append',
                            help='The column to apply the algorithm to. Defaults to all columns.')
    args = argparser.parse_args()

    # Defaults
    input_path = vars(args)['input_path'][0]
    output_path = vars(args)['output_path']
    columns = vars(args)['c'][0]
    param_m = vars(args)['m']
    param_dis = vars(args)['d']

    if len(columns) == 0:
        columns = None

    if isinstance(columns, str):
        columns = [columns]

    if not os.path.exists(input_path):
        logger.error('Input data file does not exist')
        exit()
    else:
        logger.info("Input data file: " + input_path)

    logger.info("Output file: " + output_path)

    # Load the dataset
    input_data = pd.read_csv(input_path)

    # Apply the algorithm
    output_data = suda(dataframe=input_data, max_msu=param_m, dis=param_dis, columns=columns)

    # Write the output
    output_data.to_csv(output_path, encoding='UTF-8', index=False)


if __name__ == "__main__":
    main()
