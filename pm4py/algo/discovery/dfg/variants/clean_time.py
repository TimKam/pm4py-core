'''
    This file is part of PM4Py (More Info: https://pm4py.fit.fraunhofer.de).

    PM4Py is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    PM4Py is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with PM4Py.  If not, see <https://www.gnu.org/licenses/>.
'''

from enum import Enum

import pandas as pd
from pandas.core.tools.datetimes import to_datetime

from pm4py.objects.dfg.obj import DFG
from pm4py.util import constants, exec_utils
from pm4py.util import xes_constants as xes_util
from sympy import divisors

class Parameters(Enum):
    ACTIVITY_KEY = constants.PARAMETER_CONSTANT_ACTIVITY_KEY
    CASE_ID_KEY = constants.PARAMETER_CONSTANT_CASEID_KEY
    TIMESTAMP_KEY = constants.PARAMETER_CONSTANT_TIMESTAMP_KEY


CONST_AUX_ACT = 'aux_act_'
CONST_AUX_CASE = 'aux_case_'
CONST_COUNT = 'count_'


from pm4py.discovery import discover_dfg_typed


def apply(log: pd.DataFrame, parameters=None, 
          reduction_to_abs=None, reduction_to_perc=None
          ):
    if parameters is None:
        parameters = {}

    aggregation = parameters['aggregation'] if 'aggregation' in parameters.keys() else 'avg'
    loop_handling = parameters['loop_handling'] if 'loop_handling' in parameters.keys() else 'avg'
    round_to = parameters['round_to'] if 'round_to' in parameters.keys() else 's'

    act_key = exec_utils.get_param_value(
        Parameters.ACTIVITY_KEY, parameters, xes_util.DEFAULT_NAME_KEY)
    cid_key = exec_utils.get_param_value(
        Parameters.CASE_ID_KEY, parameters, constants.CASE_ATTRIBUTE_GLUE)
    time_key = exec_utils.get_param_value(
        Parameters.TIMESTAMP_KEY, parameters, xes_util.DEFAULT_TIMESTAMP_KEY)

    '''sort the values according to cid and timekey'''
    df = log.sort_values([cid_key, time_key]).loc[:, [cid_key, act_key, time_key]].reset_index()

    '''mapping each case ID to the first time stamp of that case.'''
    grouped = df.groupby(cid_key)
    all_cases = df[cid_key].unique()

    time_dictionary_list = {}

    for case in all_cases:
        current_group = grouped.get_group(case)
        current_group = current_group.sort_values([cid_key, time_key]).loc[:, [cid_key, act_key, time_key]].reset_index()
        all_act = current_group[act_key].unique()
        init_timestamp = current_group[time_key][0]
        '''deal with loops in a case'''
        for act in all_act:
            df1 = current_group[current_group[act_key] == act]
            match loop_handling:
                case 'avg':
                    agg_time = df1[time_key].mean()
                case 'median':
                    agg_time = df1[time_key].median()
                case 'max':
                    agg_time = df1[time_key].max()
                case 'min':
                    agg_time = df1[time_key].min()
            rel_time = agg_time - init_timestamp if agg_time - init_timestamp >= pd.Timedelta(0, 's') else pd.Timedelta(0, 's')
            if act not in time_dictionary_list.keys():
                time_dictionary_list[act] = []
                time_dictionary_list[act].append(rel_time)
            else:
                time_dictionary_list[act].append(rel_time)

    keys = time_dictionary_list.keys()
    time_dictionary = {}
    for key in keys:
        match aggregation:
            case 'avg':
                agg = pd.to_timedelta(pd.Series(time_dictionary_list[key])).mean()
            case 'median':
                agg = pd.to_timedelta(pd.Series(time_dictionary_list[key])).median()
            case 'max':
                agg = pd.to_timedelta(pd.Series(time_dictionary_list[key])).max()
            case 'min':
                agg = pd.to_timedelta(pd.Series(time_dictionary_list[key])).min()
        
        time_dictionary[key] = agg.round(round_to)
    
    # reduce the number of nodes in the timeline
    try:
        if reduction_to_abs is not None:
            if reduction_to_perc is not None:
                print("Message: 'reduction_to_abs' overrides 'reduction_to_perc'.")
            time_dictionary = reduce_number_of_nodes_in_timeline(time_dictionary, reduction_to_abs=reduction_to_abs).copy()
        elif reduction_to_perc is not None: 
            time_dictionary = reduce_number_of_nodes_in_timeline(time_dictionary, reduction_to_perc=reduction_to_perc).copy()
    except Exception as e:
        print(e)

    return time_dictionary

def find_closest_natural_divisor_to_target_number(
        start_number:int, target_number=None, target_percentage=None) -> int:
    """
    Finds the divisor in a list of divisors of natural numbers 
    that is closest to a target number (absolute or percentage).
    
    Args:
        start_number (int): start number
        target_number (int): target number (absolute)
        target_percentage (float): target number (percentage)
    
    Returns:
        closest_divisor (integer): result of the closest divisor
    """
    divisor_list = divisors(start_number)
    if target_percentage:
        if target_percentage >= 0 and target_percentage <= 1:
            target_number = start_number*target_percentage
    closest_divisor = min(divisor_list, key=lambda divisor: abs(divisor - target_number))
    return closest_divisor

def replace_subsequent_sorted_values_in_dictionary(dictionary:dict, next_values:int) -> dict:
    """
    Replaces the values of a dictionary with the same value 
    for the next 'next_values' values starting with the highest.
    
    Args:
        dictionary (dict): Dictionary to be updated
        next_values (int): Number of values
    Returns:
        dictionary_update (dict): The updated dictionary
    """
    dictionary_update = dict(sorted(dictionary.items(), key=lambda item: item[1], reverse=True)).copy()
    tracker = 0
    value_replacement = 0
    for key in dictionary_update.keys():
        if tracker == 0:
            tracker = next_values + 0
            value_replacement = dictionary_update[key]
        dictionary_update[key] = value_replacement
        tracker -= 1
    return dictionary_update

def reduce_number_of_nodes_in_timeline(dictionary:dict, reduction_to_abs=None, 
                                       reduction_to_perc=None) -> dict:
    """
    Reduces the number of nodes in the timeline dictionary with time folding 
    by replacing subsequent values starting from the highest.
    
    Args:
        dictionary (dict): Timeline dictionary (must be sorted)
        reduction_to_abs (int): Reduction to an absolute number of nodes
        reduction_to_perc (float): Reduction to a percentage of the number of total nodes
    
    Returns:
        dfg_time_update (dict): Updated timeline dictionary
    """
    timenodes = len(dictionary)
    num_nodes = find_closest_natural_divisor_to_target_number(
        timenodes, target_number=reduction_to_abs, target_percentage=reduction_to_perc)
    time_replacement_num = int(timenodes/num_nodes)
    dfg_time_update = replace_subsequent_sorted_values_in_dictionary(dictionary, time_replacement_num)
    return dfg_time_update