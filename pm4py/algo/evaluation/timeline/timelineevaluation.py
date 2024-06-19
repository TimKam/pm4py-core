# packages
import pm4py
import pandas as pd
import graphviz
import math
import numpy as np
import os
from pm4py.algo.discovery.dfg.variants import clean_time
from pm4py.algo.discovery.dfg.variants import timelinelogprocessing as dft
from pm4py.visualization.dfg import visualizer as dfg_visualizer
from pm4py.visualization.dfg.variants import timeline as timeline_gviz_generator
from sklearn.linear_model import LinearRegression 
from scipy.stats import pearsonr, spearmanr

####EVALUATION SCRIPT
###

def evaluationScript(
    folder_path_import='evaluation/evaluation_data/import/', 
    folder_path_export='evaluation/evaluation_data/export/', 
    filter_cases = 0, filter_variants_k = 0, filter_variants_per = 0):
    """ Evaluates multiple logs and exports the outcome.
    """
    log_files = importMultipleEventLogs(folder_path=folder_path_import)
    print('== Evaluation starts ==')
    eva_num = 1
    for log in log_files:
        print('\n-----------------')
        print(f'start evaluation: {eva_num}')
        print('-----------------')
        print('Event log: ', log)
        print('- Import log')
        path_log = folder_path_import+log
        df = pm4py.read_xes(path_log)
        num_variants = len(pm4py.get_variants(df))
        path_statistics = folder_path_export+'stats_'+log[:-4]+'_'+str(filter_variants_per)+'_'+str(num_variants)+'.xlsx'
        try:
            
            df = dft.simplifyLog(df, 
                                 lifecycle_activities = True, 
                                 filter_cases = filter_cases, 
                                 filter_variants_k = filter_variants_k,
                                 filter_variants_per = filter_variants_per)
            print('- Import succesful')
            print('- Evaluate log')
            statistics_table = timelineEvaluationScript(df)
            print('- Evaluate succesful')
            print('- Export statistics')
            statistics_table.to_excel(path_statistics)
            print('- Export succesful')
        except:
            print('Problem occured: ', log)
        print('-----------------')
        print(f'end evaluation: {eva_num}')
        print('-----------------')
        eva_num += 1
    print('\n== Evaluation ends == ')
    return None

def importMultipleEventLogs(folder_path='evaluation/evaluation_data/import'):
    """ Returns a list of .xes file names in a folder.
    """
    files = os.listdir(folder_path)
    log_files = []
    for f in files:
        if f.endswith(".xes"):
            log_files.append(f)
    return log_files

def timelineEvaluationScript(df: pd.DataFrame, dataframe=True):
    """ Evaluation script for timeline-based discovery.
    Input is a log. It is discovered in different versions, then evaluated. 
    
    Parameters
    ==========
    df: Event log 
    dataframe: True (default) returns dataframe, False returns dictionary
    """
    logname = 'log'
    # various preprocessing strategies of the same log
    DF_dict = multipleLogPreProcessing(logname = logname, df = df)
    
    # evaluation
    STAT_dict = multipleLogEvaluation(logname = logname, DF_dict = DF_dict)
    
    if dataframe:
        STAT_dict = pd.DataFrame.from_dict(STAT_dict, orient='index')
    
    return STAT_dict

# Note: Add new pre-processing algorithms here.
def multipleLogPreProcessing(logname: str, df: pd.DataFrame):
    """ Returns a dictionary of pre-processed logs.
    """
    df0_logname = logname + '_original'
    df0 = df.copy()
    print('-- original log created')
    df1_logname = logname + '_dafsa'
    df1 = dft.dafsaDFG(df.copy())
    print('-- dafsa log created')
    df2_logname = logname + '_timerule'
    df2 = dft.timeruleDFG(df.copy())
    print('-- timerule log created')
    df3_logname = logname + '_hybrid'
    df3 = dft.timeruleDFG(df1.copy())
    print('-- hybrid log created')
    DF_dict = {
        df0_logname: df0, df1_logname: df1, 
        df2_logname: df2, df3_logname: df3}
    return DF_dict

def multipleLogEvaluation(logname: str, DF_dict: dict):
    """ Evaluates a set of logs in a dictionary
    and returns a table of statistics. 
    """
    statistics_table = {}
    df_original = logname + '_original'
    
    for log in DF_dict.keys():
        # factors
        num_a = countActivities(DF_dict[df_original]) # amount of activities before changes
        num_v = countSeqVariants(DF_dict[log])
        recurrency_dict = countRecurrentActivities(DF_dict[log])

        # metrics (based on new log)
        dfg_size = calculateDFGSize(DF_dict[log])
        contra_num, dfs_all_num = calculateDFGContradictions(DF_dict[log], selfloops=True)
        contra_num_noloops, dfs_all_num_noloops = calculateDFGContradictions(DF_dict[log], selfloops=False)
        #data_s, b0_s, b1_s, r2_s, coordination_layout_s = calculateDFGLayoutCorrelations(DF_dict[log], dfg_type='standard')
        data_t, coordination_layout_t, r_pearcorr_xy, p_pearcorr_xy, r_spearcorr_xy, p_spearcorr_xy = calculateDFGLayoutCorrelations(DF_dict[log], dfg_type='timeline', distance_type='xy_euclidian')
        data_t, coordination_layout_t, r_pearcorr_y, p_pearcorr_y, r_spearcorr_y, p_spearcorr_y = calculateDFGLayoutCorrelations(DF_dict[log], dfg_type='timeline', distance_type='y_euclidian')
        data_t, coordination_layout_t, r_pearcorr_y0, p_pearcorr_y0, r_spearcorr_y0, p_spearcorr_y0 = calculateDFGLayoutCorrelations(DF_dict[log], dfg_type='timeline', distance_type='y')
        
        # create a dictionary of all statistics
        statdict = createStatDict(logname=log, 
                                  amount_activities=num_a, 
                                  amount_seqvariants=num_v,
                                  recurr_events=recurrency_dict['recurr_events'], 
                                  recurr_sequences=recurrency_dict['recurr_sequences'], 
                                  number_events=recurrency_dict['number_events'], 
                                  number_sequences=recurrency_dict['number_sequences'], 
                                  share_events=recurrency_dict['share_events'], 
                                  share_sequences=recurrency_dict['share_sequences'],
                                  dfgsize_nodes = dfg_size[0],
                                  dfgsize_edges = dfg_size[1],
                                  #df_relations = dfs_all_num,
                                  contradictions_includesselfloops = contra_num,
                                  contradictions_excludesselfloops = contra_num_noloops, 
                                  #dfg_standard_b0 = b0_s, 
                                  #dfg_standard_b1 = b1_s, 
                                  #dfg_standard_r2 = r2_s, 
                                  
                                  # XY-axis euclidian
                                  r_pearcorr_xy = r_pearcorr_xy, 
                                  p_pearcorr_xy = p_pearcorr_xy,
                                  r_spearcorr_xy = r_spearcorr_xy, 
                                  p_spearcorr_xy = p_spearcorr_xy,
                                  
                                  # Y-axis euclidian
                                  r_pearcorr_y = r_pearcorr_y, 
                                  p_pearcorr_y = p_pearcorr_y,
                                  r_spearcorr_y = r_spearcorr_y, 
                                  p_spearcorr_y = p_spearcorr_y, 
                                  
                                  # Y-axis simple distance
                                  r_pearcorr_y0 = r_pearcorr_y0, 
                                  p_pearcorr_y0 = p_pearcorr_y0,
                                  r_spearcorr_y0 = r_spearcorr_y0, 
                                  p_spearcorr_y0 = p_spearcorr_y0)
        statistics_table.update(statdict)
        print(f'-- log {log} evaluated')
    return statistics_table

def createStatDict(logname: str, **stats):
    """ Returns a dictionary (row) with statistics from an evaluation.
    """
    statdict = {logname: stats}
    return statdict

####FACTORS
###
def countActivities(df):
    '''Return number of activity types in a log.
    '''
    ACTIVITY_COL = 'concept:name'
    activity_list = df[ACTIVITY_COL].unique()
    num_a = len(activity_list)
    #print(f"no. of activities: {num_a}")
    return num_a

def countSeqVariants(df):
    '''Return number of sequence variants in a log.
    '''
    VARIANT_COL = 'variant:concept:name'
    num_v = len(dft.variantSequences(df)[VARIANT_COL].unique())
    #print(f"no. of seq. variants: {num_v}")
    return num_v

def countRecurrentActivities(df):
    '''Return number of recurrencies of the activities in a log per case.
    '''
    ACTIVITY_COL = 'concept:name'
    CASE_COL = 'case:concept:name'
    TIME_COL = 'time:timestamp'
    #activity_list = df[ACTIVITY_COL].unique()

    #new
    df_group_reccurency = df[df.groupby(
        [CASE_COL, ACTIVITY_COL])[ACTIVITY_COL].transform('count') > 1]
    #recurr_events = len(df_group_reccurency) # number of all recurrent events
    recurr_events = len(df_group_reccurency.groupby(ACTIVITY_COL)) # number of event types with at least one recurrency in some sequence
    recurr_sequences = len(df_group_reccurency.groupby([CASE_COL])) # number of sequences with at least one recurrency
    #number_events = len(df[ACTIVITY_COL])
    number_events = len(df[ACTIVITY_COL].unique())
    number_sequences = len(df[CASE_COL].unique())
    share_events = round(recurr_events / number_events, 3)
    share_sequences = round(recurr_sequences / number_sequences, 3)
    
    recurrency_dict = {
        'recurr_events': recurr_events, 'recurr_sequences': recurr_sequences, 
        'number_events': number_events, 'number_sequences': number_sequences, 
        'share_events': share_events, 'share_sequences': share_sequences}
    
    #print(f'share of events:{recurr_events} / {number_events} = {share_events}')
    #print(f'share of sequences:{recurr_sequences} / {number_sequences} = {share_sequences}')
    
    #for activity in activity_list:
    #    df_group = df[df[ACTIVITY_COL]==activity].groupby(CASE_COL)[TIME_COL].count() > 1
    #    recurrent_activities += sum(df_group)
    
    
    
    # just recurrencies
    #for activity in activity_list:
    #    df_group = df[df[ACTIVITY_COL]==activity].groupby(CASE_COL)[TIME_COL].count()
    #    recurrency_dict[activity] = dict(df_group.aggregate([
    #        'min', 'max', 'mean', 'median']))
    #return recurrency_dict
    return recurrency_dict

####METRICS
###

#TODO: Continue with this metric
def calculateLayoutCorrelations():
    '''Return the correlations between arc length and time differences in a DFG.
    '''
    return None

def calculateDFGSize(df):
    '''Return the number of nodes and edges of a DFG. 
    Does not count the start and end activities.
    '''
    nodes = countActivities(df)
    edges = len(pm4py.discover_directly_follows_graph(df)[0])  
    dfg_size = (nodes, edges)
    #print(f"no. of nodes: {nodes}")
    #print(f"no. of edges: {edges}")
    return dfg_size

def calculateDFGContradictions(df, selfloops=True):
    '''Return number of contradictions to timeline and 
    directly follows relations between directly follows relations.
    selfloops = 'true' includes selfloops as contradictions (default)
    '''
    timeline_comparison = dft.createDFallnext(dft.TimeLine(df).timeDict())
    #print(timeline_comparison)
    #timeline_dict = clean_time.apply(df)
    #timeline_comparison = dft.createDFallnext(timeline_dict)
    dfg, start_activities, end_activities = pm4py.discover_directly_follows_graph(df)
    dfs_all_num = len(dfg.keys())
    contra_num = 0
    for relation in dfg.keys():
        #print(relation[0])
        #print("->",relation[1])
        #print(timeline_comparison[relation[0]])
        if not {relation[1]}.issubset(timeline_comparison[relation[0]]):
            #print(True)
            # don't count selfloops if applied
            if (selfloops==False and relation[1] == relation[0]):
                None
            else:
                contra_num += 1
        else:
            #print(False)
            None
    #print("contradictions: ", contra_num)
    #print("all dfs: ", dfs_all_num)
    return (contra_num, dfs_all_num)

def calculateDFGLayoutCorrelations(
    df: pd.core.frame.DataFrame, dfg_type='standard', distance_type='xy_euclidian'
) -> tuple[dict, np.float64, np.float64, np.float64, dict]: 
    '''Calculate the linear regression and R^2 for the layout from an event log.
    
    Parameters
    ==========
    df: pandas dataframe
    dfg_type: 'standard' or 'timeline'
    distance_type: 'xy_euclidian' (for 2D), 'x_euclidian', 'x' (difference), 'y_euclidian', 'y' (difference)
    '''
    if dfg_type=='standard':
        dfg, start_activities, end_activities = pm4py.discover_directly_follows_graph(df)
        gviz = dfg_visualizer.apply(dfg)
    elif dfg_type=='timeline':
        dfg, start_activities, end_activities = pm4py.discover_dfg_typed(df)
        dfg_time = clean_time.apply(df)
        gviz = timeline_gviz_generator.apply(
            dfg, dfg_time, parameters={
                "format": "png", 
                "start_activities": start_activities,"end_activities": end_activities})
    else:
        print('Please choose DFG \'standard\' or \'timeline\'')
        return None
    
    d_index = 0 # index for distance type
    t_index = 5 # index for time difference
    if distance_type == 'xy_euclidian':
        None
    elif distance_type == 'x_euclidian':
        d_index = 1
    elif distance_type == 'y_euclidian':
        d_index = 2
    elif distance_type == 'x':
        d_index = 3
    elif distance_type == 'y':
        d_index = 4
    else:
        print('You have to choose distance_type: xy_euclidian, x_euclidian, y_euclidian, x, or y.')
    
    # extract coordinates from gviz
    #print(type(df['concept:name']))
    activities_list = df['concept:name'].unique()
    coordination_layout, df_layout = extractGvizCoordinates(gviz, activities_list)
    
    # calculate distances based on coordinates and timeline
    #timeline_dict = dft.TimeLine(df).timeDict()
    timeline_dict = clean_time.apply(df) # time nodes on timeline-based approach
    coordination_layout = calculateNodeDistances(coordination_layout, timeline_dict)
    
    # calculate pearson correlation
    data = {'distance': [], 'time_difference_second': []}
    for edge in coordination_layout['distances']:
        data['distance'].append(coordination_layout['distances'][edge][d_index])
        data['time_difference_second'].append(coordination_layout['distances'][edge][t_index])
    X = np.array(data['distance'])#.reshape(-1, 1)
    Y = np.array(data['time_difference_second'])#.reshape(-1, 1)
    r_pearcorr, p_pearcorr = pearsonr(X, Y)
    r_spearcorr, p_spearcorr = pearsonr(X, Y)
    #print('r_corr pear:', r_pearcorr)
    #print('r_corr spear:', r_spearcorr)
    #print('p_value pear:', p_pearcorr)
    #print('p_value spear:', p_spearcorr)
    
    # calculate a simple linear regression model (could be deleted)
    #data = {'euclidian_distance': [], 'time_difference_second': []}
    #for edge in coordination_layout['distances']:
    #    data['euclidian_distance'].append(coordination_layout['distances'][edge][0])
    #    data['time_difference_second'].append(coordination_layout['distances'][edge][1])
    #X = np.array(data['euclidian_distance']).reshape(-1, 1)
    #y = np.array(data['time_difference_second']).reshape(-1, 1)
    #lr = LinearRegression().fit(X, y)
    #b0 = lr.intercept_
    #b1 = lr.coef_
    #r2 = lr.score(X, y)
    #print(f'intercept: {b0}, coefficiant: {b1[0]}, R^2: {r2}')
    #print(data)
    #return data, b0[0], b1[0][0], r2, coordination_layout
    return data, coordination_layout, r_pearcorr, p_pearcorr, r_spearcorr, p_spearcorr

def extractGvizCoordinates(
    gviz: graphviz.graphs.Digraph, activities_list = [list]
) -> tuple[dict, dict, pd.core.frame.DataFrame]:
    """ Returns a dictionary and a pandas dataframe with coordinates of nodes from an graphviz object.
    ===== 
    Parameters
    =====
    gviz: graphiviz object 
    """
    dot_file = gviz.pipe(format='plain', encoding='utf-8')
    img_text = dot_file.split('\n')
    """
    How the file from the plain text format in graphviz is structured:
    - row1: graph scale width height
    - row2: node name x y width height label style shape color fillcolor
    - row3: edge tail head n x₁ y₁ .. xₙ yₙ [label xl yl] style color
    # stop    
    reference: https://www.graphviz.org/docs/outputs/plain/
    """

    i=0
    coordination_layout = {
        'nodes': [], 'edges': [],
        'node_names': {}, 'node_ids': {},
        'coordinates': {}}
    excluded_nodes = []
    for i in range(0, len(img_text)):
        row_info = img_text[i].split(' ')
        # Collect all node data
        if 'node' in row_info[0]:
            #print(set(row_info))
            row_type = 'node'
            node_id = row_info[1]
            node_pos_x = row_info[2]
            node_pos_y = row_info[3]
            node_activity = img_text[i].split('"')[1].rsplit(" ", 1)[0]
            if node_activity in activities_list:
                coordination_row = {'row_type': row_type,
                                    'node_id': node_id, 
                                    'node_pos_x': node_pos_x, 
                                    'node_pos_y': node_pos_y, 
                                    'node_activity': node_activity
                                   }
                coordination_layout['nodes'].append(coordination_row)
                coordination_layout['node_ids'][node_id] = node_activity
                coordination_layout['node_names'][node_activity] = node_id
            else:
                #print('node excluded:', node_activity)
                excluded_nodes.append(node_id)
                continue;
            
            
            try:
                coordination_layout['coordinates'][node_id] = [float(node_pos_x), float(node_pos_y)]
            except ValueError:
                coordination_layout['coordinates'][node_id] = 'NaN'
        # Collect all edge data
        elif 'edge' in row_info[0]:
            row_type = 'edge'
            from_node = row_info[1]
            to_node = row_info[2]
            if not set(excluded_nodes).intersection(set(row_info)):
                coordination_row = {'row_type': row_type,
                                'from_node': from_node, 
                                'to_node': to_node, 
                               }
                coordination_layout['edges'].append(coordination_row)
            else:
                None
                #print('edge excluded:', (from_node, to_node))
                continue;
                
    complete_table = coordination_layout['nodes'] + coordination_layout['edges']
    df_layout = pd.DataFrame(complete_table)
    return coordination_layout, df_layout

def calculateNodeDistances(coordination_layout: dict, timeline_dict: dict):
    """ Extends the coordination_layout with distance between edges.
    Calculates both euclidian distance and distance between timenodes for the same edge.
    """
    # calculate edge distances 
    edge_distance_eucl_dict = {}
    edge_distance_eucl_x_dict = {} # NEW
    edge_distance_eucl_y_dict = {} # NEW
    edge_distance_x_dict = {} # NEW
    edge_distance_y_dict = {} # NEW
    edge_distance_time_dict = {}
    edge_distances = {}
    for edge in coordination_layout['edges']:
        # calculate euclidian distance between two points p and q
        p = coordination_layout['coordinates'][edge['from_node']]
        q = coordination_layout['coordinates'][edge['to_node']]
        d2D = math.dist(p, q)
        d1Dx = math.dist([p[0]], [q[0]])
        d1Dy = math.dist([p[1]], [q[1]])
        edge_distance_eucl_dict[(edge['from_node'], edge['to_node'])] = d2D
        edge_distance_eucl_x_dict[(edge['from_node'], edge['to_node'])] = d1Dx
        edge_distance_eucl_y_dict[(edge['from_node'], edge['to_node'])] = d1Dy
        
        # calculate the distance between two points p and q
        d0x = p[0] - q[0]
        d0y = p[1] - q[1]
        edge_distance_x_dict[(edge['from_node'], edge['to_node'])] = d0x
        edge_distance_y_dict[(edge['from_node'], edge['to_node'])] = d0y
        
        # calculate time difference between two timenodes
        t1 = timeline_dict[coordination_layout['node_ids'][edge['from_node']]]
        t2 = timeline_dict[coordination_layout['node_ids'][edge['to_node']]]
        t_diff = (t2-t1).total_seconds()
        edge_distance_time_dict[(edge['from_node'], edge['to_node'])] = t_diff
        
        # combine all distances in one dictionary
        edge_distances[(edge['from_node'], edge['to_node'])] = (d2D, d1Dx, d1Dy, d0x, d0y, t_diff)
        
    complete_table = coordination_layout['nodes'] + coordination_layout['edges']
    coordination_layout['distances_euclidian'] = edge_distance_eucl_dict
    coordination_layout['distances_euclidian_x'] = edge_distance_eucl_x_dict
    coordination_layout['distances_euclidian_y'] = edge_distance_eucl_y_dict
    coordination_layout['distances_x'] = edge_distance_x_dict
    coordination_layout['distances_y'] = edge_distance_y_dict
    coordination_layout['distances_time'] = edge_distance_time_dict
    coordination_layout['distances'] = edge_distances
    return coordination_layout