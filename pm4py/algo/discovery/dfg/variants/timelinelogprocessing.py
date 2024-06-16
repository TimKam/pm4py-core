import pm4py
import pandas as pd
import numpy as np
from dafsa import DAFSA

def helloWorld():
    var = 'helloWorld'
    return print(var)

####TIMELINE RULES
###
# MAIN FUNCTION
def timeruleDFG(df, time_calculation='mean'):
    ''' Creates a timerule-DFG by changing the activity names accordingly. 
    '''
    # keep only necessary attributes
    #df = simplifyLog(df) 
    # create a column with relative timestamps
    df = relativeTimestamps(df)

    # 1) Divide log in sequence variants
    #print('1) Log divided')
    df = variantSequences(df)
    variants = df['variant:concept:name'].unique()
    
    # 2) Create/ Update (updates below) timeline a timeline
    #print('2) Create timelines')
    time_0 = TimeLine(df, time_calculation=time_calculation) 
    timedict_0 = time_0.timeDict() # could be implemented in the class
    time_DF_0 = createDFall(timedict_0) # could be implemented in the class
    activity_versions = time_0.activityVersions(timedict_0)

    
    # 3) Check each variant 
    for v in variants:
        #print('\n====check variant: ', v)
        # 3.a) Create a variant timeline.  
        df_variant = df[df['variant:concept:name']==v]
        time_1 = TimeLine(df_variant, time_calculation=time_calculation)
        timedict_1 = time_1.timeDict()
        time_DF_1 = createDFnext(timedict_1)
        
        # 3.b) Compare base-timeline with variant-timeline for 
        #print('-> timelineComparison()')
        activity_changes = timelineComparison(timedict_0, time_DF_0, timedict_1, time_DF_1, activity_versions)

        # Rename activities that do not adhere to timeline
        if activity_changes:
            for i in range(0, len(activity_changes)):
                old_activity = activity_changes[i][0]
                new_activity = activity_changes[i][1]
                #new code
                indices = df['concept:name'].loc[
                    (df['variant:concept:name']==v) & (df['concept:name']==old_activity)].index.to_list()
                df.loc[indices, 'concept:name'] = new_activity
                # old code
                #df['concept:name'].loc[df['variant:concept:name']==v] = df['concept:name'].loc[df['variant:concept:name']==v].replace(old_activity, new_activity)
                
        # (2) Update timeline
        time_0 = TimeLine(df, time_calculation=time_calculation) 
        timedict_0 = time_0.timeDict()
        time_DF_0 = createDFall(timedict_0)
        activity_versions = time_0.activityVersions(timedict_0)
        #print('activity_versions', activity_versions)
    return df

class TimeLine:
    '''Create and maintain timelines using a log (dataframe).
    '''
    ACTIVITY_COL = 'concept:name'
    TIMEREL_COL = 'time:relative'
    TIMELINE_DICT = {}
    ACTIVITYVER_DICT = {}
    
    def __init__(self, dataframe, time_calculation='mean'):
        self.timecalculation = time_calculation
        self.df = relativeTimestamps(dataframe)
    
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
def timelineComparison2(timedict_0, time_DF_0, timedict_1, time_DF_1, activity_versions):
    '''Compares a two timelines to find discrepancies.
    '''
    activity_changes = []
    a_change_check = ('', '')
    #print('======')
    #print('timedict_0', timedict_0)
    #print('\ntimedict_1', timedict_1)
    #print('\ntime_DF_0', time_DF_0)
    #print('\ntime_DF_1', time_DF_1)
    #print('\nactivity_versions', activity_versions)
    #print('======\n')
    # Check each relation a and b, where a < b in timeline 1.
    for activity_a in time_DF_1:
        activity_b = next(iter(time_DF_1[activity_a]))
        # needed for activity version checking
        level_skip = False
        sep = '_'
        activity_a_original = activity_a.split(sep, 1)[0]
        activity_b_original = activity_b.split(sep, 1)[0]
        activity_a_dummy = ''
        activity_b_dummy = ''
        #print('activity_a', activity_a)
        #print('activity_b', activity_b)
        
            
        # CHECK 1) a < b => A < B | a=A, b=B, a≠a', a≠a*
        if (({activity_b}.issubset(time_DF_0[activity_a]) and 
             a_change_check[0] == '') or 
            (activity_b == 'final')):
            #print('- level 1')
            level_skip = True;
        
        # CHECK 2) a' < b => A' < B | a'=A', b=B
        if (level_skip == False):
            #print('- level 2')
            # - version already exists (a=a')
            if a_change_check[0] == 'version':
                if ({activity_b}.issubset(time_DF_0[a_change_check[1]])):
                    a_change_check = ('', '')
                    level_skip = True;
            # - version must be looked up (a≠a')
            elif a_change_check[0] == '':
                for activity_a_version in activity_versions[activity_a_original]:
                    if ({activity_b}.issubset(time_DF_0[activity_a_version])):
                        activity_a_dummy = activity_a_version
                        level_skip = True;
                        break;
        
        # CHECK 3) a < b' => A < B' | a=A, b'=B', a'≠a*
        if (level_skip == False):
            #print('- level 3')
            for activity_b_version in activity_versions[activity_b_original]:
                if ({activity_b_version}.issubset(time_DF_0[activity_a])):
                    activity_b_dummy = activity_b_version
                    a_change_check = ('version', activity_b_dummy)
                    level_skip = True
                    break;
        
        # CHECK 4) a* < b => t(a*) < t(B) | b=B
        if (level_skip == False):
            #print('- level 4')
            # create a new dummy activity for a -> a* if not new
            if a_change_check[0] == 'new':
                a_change_check = ('', '')
            elif a_change_check[0] == '':
                number = len(activity_versions[activity_a_original])
                activity_a_dummy = f"{activity_a_original}_v{number}"
            
            if timedict_1[activity_a] <= timedict_0[activity_b]:
                level_skip = True
                
        # CHECK 5) a* < b' => t(a*) < t(B') | b=B
        if (level_skip == False):
            #print('- level 5')
            for activity_b_version in activity_versions[activity_b_original]:
                if (timedict_1[activity_a] <= timedict_0[activity_b_version]):
                    activity_b_dummy = activity_b_version
                    a_change_check = ('version', activity_b_dummy)
                    level_skip = True
                    break;
        
        # CHECK 6) a* < b' => t(a*) < t(B') | b=B
        if (level_skip == False):
            #print('- level 6')
            # create a new dummy activity for b -> a*
            number = len(activity_versions[activity_b_original])
            activity_b_dummy = f"{activity_b_original}_v{number}"
            a_change_check = ('new', activity_b_dummy)
        
          
        # FINALIZE: change activities
        if activity_a_dummy != '':
            #print('=> activity_a_dummy', activity_a_dummy)
            activity_changes.append((activity_a, activity_a_dummy))
            activity_a_dummy = ''
        if activity_b_dummy != '':
            #print('=> activity_b_dummy', activity_b_dummy)
            activity_changes.append((activity_b, activity_b_dummy))
            activity_b_dummy = ''
        #print('=> activity_changes', activity_changes)
            
    return activity_changes

def timelineCompariso3(timedict_0, time_DF_0, timedict_1, time_DF_1, activity_versions):
    '''Compares a two timelines to find discrepancies.
    '''
    activity_changes = []
    timediff = {'start'}
    #print(timedict_0)
    for activity in time_DF_1:
        next_activity = next(iter(time_DF_1[activity]))
          
        if (
            (time_DF_1[activity].issubset(time_DF_0[activity]) and not
            timediff.intersection(time_DF_0[next_activity]))
        ):
            None
        elif (next_activity == 'final' and not
              timediff.intersection(time_DF_0[activity])
             ):
            None
        else: 
            dummy_activity = ""
            sep = '_'
            current_activity = activity.split(sep, 1)[0]
            for version in activity_versions[current_activity]:
                if (
                    time_DF_1[activity].issubset(time_DF_0[version]) and not
                    timediff.intersection(time_DF_0[version])
                   ):
                    dummy_activity = version
                    break;
            if dummy_activity == "":
                number = len(activity_versions[current_activity])
                dummy_activity = f"{activity}_v{number}"
            activity_changes.append((activity, dummy_activity))
        
        if next_activity != 'final':
            timediff.update(time_DF_0[activity].difference({next_activity}))
            timediff.update({activity})
    return activity_changes

def timelineComparison(timedict_0, time_DF_0, timedict_1, time_DF_1, activity_versions):
    '''Compares a two timelines to find discrepancies.
    '''
    activity_changes = []
    previous_activity = ''
    for activity in time_DF_1:
        if (
            # Check 1: Leave start activity as is
            previous_activity == ''
        ):
            None
        elif (
            'final' in time_DF_1[activity] and
            timedict_1[previous_activity]<=timedict_0[activity]
            
        ):
            None
        elif (
            time_DF_1[activity].issubset(time_DF_0[activity]) and not
            {previous_activity}.issubset(time_DF_0[activity])
        ):
            None
        else: 
            dummy_activity = ""
            sep = '_'
            current_activity = activity.split(sep, 1)[0]
            for version in activity_versions[current_activity]:
                if (
                    time_DF_1[activity].issubset(time_DF_0[version]) and not 
                    {previous_activity}.issubset(time_DF_0[version])
                   ):
                    dummy_activity = version
                    break;
            if dummy_activity == "":
                number = len(activity_versions[current_activity])
                dummy_activity = f"{activity}_v{number}"
            activity_changes.append((activity, dummy_activity))
        previous_activity = activity
    return activity_changes

# ToDo delete and use createDFallnext() instead?
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

def createDFallnext(timedict: dict):
    '''Returns a dictionary with the activities that follow one activity in timeline dictionary.
    Warning! Dictionary must be sorter according to time/occurency.
    '''
    activity_list = list(timedict.keys())
    timeline_relations = {}
    for i in range(0, len(activity_list)):
        if i+1 < len(activity_list):
            timeline_relations[activity_list[i]] = set(activity_list[i+1:len(activity_list)])
        else:
            timeline_relations[activity_list[i]] = {'final'}
    return timeline_relations


#### DAFSA
###
# MAIN FUNCTION
def dafsaDFG(df):
    ''' Creates a DAFSA-DFG by changing the activity names accordingly. 
    '''
    # 1) Rename activities AND 2) Create sequence strings
    sequences, a_seq_dict = sequenceCreator(df)
    #print(sequences)
    #print('\n', a_seq_dict)
    # 3) Create a DAFSA graph
    d = DAFSA(sequences) # function sorts the sequences
    
    # 4) Find paths in graph and create node id sequences for each case
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
    
    # 5) Create dummy activities 
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
    '''Give each activity in a log a unique character.
    '''
    # strings of possible characters
    char_upper_case = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    char_lower_case = 'abcdefghijklmnopqrstuvwxyz'
    char_num = '0123456789'
    char_sym = '~!@#$%^&*()+={[}]|\:;<,>.?/' # excluded: ' " _ ` -
    # combine all strings into one
    charid = char_upper_case+char_lower_case+char_num+char_sym
    
    # add a new column with one-character id's for each activity
    activity_list = df["concept:name"].unique()
    activity_dict = {}
    if len(activity_list) > len(charid):
        print('There are too many activity types than characters.\n'
        'Please consider adding more characters or '
        'reducing the number of unique activity types.')
        return None
    a_num = 0
    for activity in activity_list:
        activity_dict[activity] = charid[a_num]
        a_num += 1
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

def simplifyLog(df, lifecycle_activities=False, filter_cases = 0):
    '''simplifies a log by keeping only necessary attributes.
    '''
    #keep following columns
    case = 'case:concept:name'
    activity = 'concept:name'
    time = 'time:timestamp'
    lifecycle = 'lifecycle:transition'
    # create new activities with transitions
    if (lifecycle_activities == True and 
        lifecycle in df.columns):
        if len(df[lifecycle].unique()) > 1:
            df[activity] = df[activity] + '-' + df[lifecycle]
        else:
            print('Message: No transition-activity were be created. Only one type of lifecycle transition in log.')
    else:
        print('Message: No transition-activity were be created. No transition column.')
    # keep only k amount cases in the log
    if filter_cases > 0:
        case_list = df[case].unique()[0:filter_cases]
        df = pm4py.filter_event_attribute_values(df, case, case_list, level="case", retain=True).copy()
        
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