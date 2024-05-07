import pm4py
import pandas as pd

from dafsa import DAFSA

####TIMELINE RULES
###
# MAIN FUNCTION
def timeruleDFG(df, time_calculation='mean'):
    ''' Creates a timerule-DFG by changing the activity names accordingly.'''
    # keep only necessary attributes
    #df = simplifyLog(df) 
    # create a column with relative timestamps
    df = relativeTimestamps(df)
    
    time_0 = TimeLine(df, time_calculation=time_calculation) 
    timedict_0 = time_0.timeDict() # could be implemented in the class
    time_DF_0 = createDFall(timedict_0) # could be implemented in the class
    activity_versions = time_0.activityVersions(timedict_0)
    # compare variant sequences with timeline
    df = variantSequences(df)
    variants = df['variant:concept:name'].unique()
    for v in variants:
        # ToDo: new function:
        #print('\nCheck variant:', v)
        #print('-----------------')
        # initialize the timeline for the variant
        df_variant = df[df['variant:concept:name']==v]
        time_1 = TimeLine(df_variant, time_calculation=time_calculation)
        timedict_1 = time_1.timeDict()
        time_DF_1 = createDFnext(timedict_1)

        activity_changes = timelineComparison(timedict_0, time_DF_0, timedict_1, time_DF_1, activity_versions)

        # rename activities that do not adhere to timeline
        if activity_changes:
            for i in range(0, len(activity_changes)):
                old_activity = activity_changes[i][0]
                new_activity = activity_changes[i][1]
                df['concept:name'].loc[df['variant:concept:name']==v] = df['concept:name'].loc[df['variant:concept:name']==v].replace(old_activity, new_activity)
        # update timeline with new activities
        time_0 = TimeLine(df, time_calculation=time_calculation) 
        timedict_0 = time_0.timeDict()
        time_DF_0 = createDFall(timedict_0)
        activity_versions = time_0.activityVersions(timedict_0)
        #print('activity_versions', activity_versions)
    return df

class TimeLine:
    '''Create and maintain timelines using a log (dataframe).'''
    ACTIVITY_COL = 'concept:name'
    TIMEREL_COL = 'time:relative'
    TIMELINE_DICT = {}
    ACTIVITYVER_DICT = {}
    
    def __init__(self, dataframe, time_calculation='mean'):
        self.timecalculation = time_calculation
        self.df = dataframe
    
    def timeDict(self):
        if self.TIMEREL_COL in self.df.columns:
            if self.timecalculation == 'mean':
                self.TIMELINE_DICT = self.df.groupby([self.ACTIVITY_COL])[self.TIMEREL_COL].mean().to_dict()
            elif self.timecalculation == 'median':
                self.TIMELINE_DICT = self.df.groupby([self.ACTIVITY_COL])[self.TIMEREL_COL].median().to_dict()
            else: 
                raise ValueError("Time calculation did not succeed. Choose \"mean\" or \"median\".")
        else: 
            raise KeyError("Dataframe does not contain relative timestamps.")
        #sort dictionary according to times
        self.TIMELINE_DICT = dict(sorted(self.TIMELINE_DICT.items(), key=lambda item: item[1]))
        return self.TIMELINE_DICT
    
    def activityVersions(self, timedict: dict):
        ''' Creates a dictionary with activity versions. Input is a dictionary of activities.
        '''
        activities = list(timedict.keys())
        activities_original = self.keepUniqueActivities(activities)
        self.ACTIVITYVER_DICT = {}
        for a in activities_original:
            self.ACTIVITYVER_DICT[a] = self.keepSingleActivities(activities, a)
        return self.ACTIVITYVER_DICT
    
    def splitString(self, word: str, sep='_'):
        '''Splits string at first occurance of a symbol.
        '''
        newword = word.split(sep, 1)[0]
        return newword

    def keepUniqueActivities(self, activities: list):
        '''Creates a list of unique activities from a list activities.
        '''
        activities_unique = []
        for a in activities:
            newactivity = self.splitString(a)
            if newactivity not in set(activities_unique):
                activities_unique.append(newactivity)
        return activities_unique

    def keepSingleActivities(self, activities: list, activity: str):
        '''Creates a list of all versions of an activity from a list activities.
        '''
        activities_unique = []
        for a in activities:
            newactivity = self.splitString(a)
            if newactivity == activity:
                activities_unique.append(a)
        activities_unique.sort()
        return activities_unique
    
def variantSequences(df):
    '''Adds sequence variants to event log (dataframe).
    '''
    df = df.sort_values('case:concept:name', ascending=True)
    variants = pm4py.get_variants(df)
    variants = list(variants.keys())
    variantseq_dict = {}
    variantnum_dict = {}
    for i in range(0, len(variants)):
        df_variant = pm4py.filter_variants(df, [variants[i]])
        for case in df_variant['case:concept:name'].unique():
            variantseq_dict[case] = variants[i]
            variantnum_dict[case] = i
    df['variant:sequence'] = df['case:concept:name'].map(variantseq_dict)
    df['variant:concept:name'] = df['case:concept:name'].map(variantnum_dict)
    return df

def createDFall(timedict: dict):
    '''Returns a dictionary with all activites that follow one activity in a sequence.
    Warning! Dictionary must be sorter according to time/occurency.
    '''
    sequence = list(timedict.keys())
    DF_dict = {}
    for activity in sequence.copy():
        del sequence[0]
        if len(sequence) != 0:
            DF_dict[activity] = set(sequence)
        else: 
            DF_dict[activity] = 'final'
    return DF_dict

# ToDo: Combine timedict and timedf. Just use the function for that.
def timelineComparison(timedict_0, time_DF_0, timedict_1, time_DF_1, activity_versions):
    '''Compares a two timelines to find discrepancies.
    '''
    
    activity_changes = []
    previous_relation = ''
    for relation in time_DF_1:
        try: 
            #print(f'\n(previous act.: {previous_relation} ({timedict_1[previous_relation]}) / ({timedict_0[previous_relation]}))')
            1+1
            None
        except: 
            None
        #print(f'Activity: {relation} ({timedict_1[relation]}) / ({timedict_0[relation]})')
        try:
            #print(f'- next: {time_DF_1[relation]} ({timedict_1[relation]})/({timedict_0[list(time_DF_1[relation])[0]]})')
            1+1
            None
        except:
            #print(f'final')
            None
        try:
            #print(f'rule DF: {time_DF_1[relation].issubset(time_DF_0[relation])}, rule Time DF: {timedict_0[previous_relation]<=timedict_0[relation]}')
            1+1
            None
        except:
            None
        #print('- ', time_DF_1[relation].issubset(time_DF_0[relation]))
        
        # check if there are any discrepancies
        if (
            #first activity
            timedict_1[relation] == pd.Timedelta('0D')
        ):
            #print("=> start activity")
            None
        elif (
            #final activity no discrepancies
            'final' in time_DF_1[relation] and
            timedict_1[previous_relation]<=timedict_0[relation] and
            timedict_0[previous_relation]<=timedict_0[relation]
        ):
            #print("=> final activity, no prob")
            None
        elif (
            timedict_0[previous_relation]<=timedict_0[relation] and
            time_DF_1[relation].issubset(time_DF_0[relation]) and not
            {previous_relation}.issubset(time_DF_0[relation])
        ):
            #print("=> OK")
            None
        else: 
            #print("----DISCREPANCY----")
            new_activity = ""
            sep = '_'
            new_relation = relation.split(sep, 1)[0]
            #print('check versions:')
            for version in activity_versions[new_relation]:
                if (
                    timedict_0[previous_relation]<=timedict_0[version] and
                    time_DF_1[relation].issubset(time_DF_0[version]) and not 
                    {previous_relation}.issubset(time_DF_0[version])
                   ):
                    #print("- previous version works")
                    new_activity = version
                    #print('new_activity (for version):', new_activity)
                    break;
            if new_activity == "":
                #print("- new version created")
                #print('activity_versions: ', activity_versions[new_relation])
                number = len(activity_versions[new_relation])
                #print('number: ', number)
                new_activity = f"{relation}_v{number}"
            #print(f'=> new activity: {new_activity}')
            activity_changes.append((relation, new_activity))
            #result = time_DF_1[relation].issubset(time_DF_0[relation])
            #activity_changes.append((relation, f'{relation}_v0'))
        previous_relation = relation
    return activity_changes

def createDFnext(timedict: dict):
    '''Returns a dictionary with the activity that follows one activity in a sequence.
    Warning! Dictionary must be sorter according to time/occurency.
    '''
    sequence = list(timedict.keys())
    DF_dict = {}
    for activity, i in zip(sequence.copy(), range(0, len(sequence))):
        try:
            DF_dict[activity] = {sequence[i+1]}
        except IndexError:
            DF_dict[activity] = {'final'}
    return DF_dict




#### DAFSA
###
# MAIN FUNCTION
def dafsa_dfg(df):
    ''' Creates a DAFSA-DFG by changing the activity names accordingly. 
    '''
    sequences, a_seq_dict = sequenceCreator(df)
    d = DAFSA(sequences)
    # create new sequences based on node ids
    a_new_seq_dict = {}
    for case in a_seq_dict:
        new_sequence = []
        # Start at the root
        node = d.lookup_nodes[0]
        # If we can follow a path, it is valid, otherwise return None
        cum_weight = 0
        for token in a_seq_dict[case]:
            if token not in node.edges:
                None
            cum_weight += node.edges[token].weight
            node = node.edges[token].node
            new_sequence.append(str(node.node_id))
        # Check if the last node is indeed a final one (so we don't
        # match prefixes alone)
        if not node.final:
            None
        a_new_seq_dict[case] = new_sequence
        
    # attach node number to index in dataframe
    index_seq_dict = {}
    for case in df['case:concept:name']:
        #print(f"case:{case}")
        seq_number=0
        for index in df.loc[df['case:concept:name']==case].index:
            #print(index)
            #print(a_new_seq_dict[case][seq_number])
            index_seq_dict[index] = a_new_seq_dict[case][seq_number]
            seq_number+=1
    # add to dataframe
    df['concept:name:nodeid'] = df.index.to_series().map(index_seq_dict)
    
    # create a new activity names based on activity name and node id
    new_col_dict = {}
    for index in df.index:
        new_name = df['concept:name'].loc[index]+"-N"+df['concept:name:nodeid'].loc[index]
        new_col_dict[index] = new_name
    # add to dataframe
    df['concept:name:namenode'] = df.index.to_series().map(new_col_dict)
    # new log
    df = df.rename(columns={
        "concept:name": "concept:name_old", 
        "concept:name:namenode": "concept:name"})
    return df

def sequenceCreator(df):
    ''' Create sequence strings for each case in a log.
    '''
    df = activityRenaming(df)
    # create "words" with case activites
    a_seq_dict = {} 
    for case in df["case:concept:name"]:
        a_sequence = df["concept:name:num"].loc[df["case:concept:name"] == case].sum()
        a_seq_dict[case] = a_sequence
    sequences = list(a_seq_dict.values())
    return sequences, a_seq_dict

def activityRenaming(df):
    '''Give each activity in a log a number.
    '''
    activity_dict = {}
    a_num = 0
    for activity in df["concept:name"].unique():
        a_num += 1
        activity_dict[activity] = str(a_num)
    df['concept:name:num'] = df["concept:name"].map(activity_dict)
    return df

####GENERAL
###
def relativeTimestamps(df):
    '''adds relative timestamps to dataframe (log).
    '''
    df = df.sort_values(by='time:timestamp')
    start_times_dict = {}
    cases = df['case:concept:name'].unique()
    for case in cases:
        start_time = df['time:timestamp'].loc[df['case:concept:name'] == case].iloc[0]
        start_times_dict[case] = start_time
        df['time:relative'] = (
            df['time:timestamp'] 
            - df['case:concept:name'].map(start_times_dict)
        )
    return df

def simplifyLog(df):
    '''simplifies a log by keeping only necessary attributes.
    '''
    #keep following columns
    case = 'case:concept:name'
    activity = 'concept:name'
    time = 'time:timestamp'
    #filter log
    keep_columns = [case, activity, time]
    df = df.filter(keep_columns)
    return df

def recurrentActivities(df):
    '''Makes every activity in a log unique 
    by adding a version number to recurrent activities.
    '''
    cases = df['case:concept:name'].unique()
    for case in cases:
        df2 = df[df['case:concept:name'] == case]
        group = df2[df2.duplicated(['concept:name'], keep=False)].groupby('concept:name').count()
        duplicates = group.index.to_list()
        duplicates.sort()
        for activity in duplicates:
            num = 0
            indices = df[(df['case:concept:name'] == case)&(df['concept:name']==activity)].index.to_list()
            for i in indices:
                if num > 0:
                    df['concept:name'].loc[i] = activity + '_' + str(num)
                num+=1
    return df