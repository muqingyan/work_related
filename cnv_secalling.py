import pandas as pd
import numpy as np
import renasight_dq_zscore_init as initcnv
from itertools import groupby
from operator import itemgetter

def recal_SE(df_normal,sample,gene):
    df_normal_check = df_normal[(df_normal['Sample'] == sample) & (df_normal['gene'] == gene)]
    m = df_normal_check.dq.mean(skipna=True)
    s = df_normal_check.dq.std(ddof=0, skipna=True)
    l = m - (3 * s)
    h = m + (3 * s)
    q1 = np.quantile(df_normal_check['dq'], 0.25)
    q3 = np.quantile(df_normal_check['dq'], 0.75)
    iqr = q3 - q1
    ub = q3 + (1.5 * iqr)
    lb = q1 - (1.5 * iqr)
    return l,h,ub,lb

def assign_se_type(row,l,h,lb,ub):
    if row['consensus_CNV'] == 'deletion' and row['dq'] <= l and row['dq'] < lb:
        return 'SE'
    elif row['consensus_CNV'] == 'duplication' and row['dq'] >= h and row['dq'] > ub:
        return 'SE'
    else:
        return np.nan
def tell_se_me(exons):
    exon_cl=[]
    for k, g in groupby(enumerate(exons), lambda ix:ix[0] -ix[1]):
        exon_cl.append(list(map(itemgetter(1),g)))
    return exon_cl

def remove_gene(df_count_adj,df_dq_filterCV):
    # df_dq = pd.read_csv(dq015path)
    df_dq = df_dq_filterCV

    print("Past CV0.15 filter original DQ",df_dq.head())
    # print("GET 10287832,'NPHP1:exon_1'",df_dq[(df_dq['Sample']==10287832)]['NPHP1:exon_1'])

    df_dq2 = pd.melt(df_dq, id_vars=['Sample'], var_name='target', value_name='dq')
    df_dq2[['gene', 'exon']] = df_dq2.target.str.split(':exon_', expand=True)
    print(df_dq2)
    print(df_count_adj)
    df_dq_adj=df_dq2.merge(df_count_adj,on=['Sample','gene'],how='outer')
    print("after merge dq and adj_cnv_type",df_dq_adj)

    # print("after merge dq and adj_cnv_type",df_dq_adj[(df_dq_adj['Sample']==12697070) & (df_dq_adj['gene']=='CHRM3')])
    # print("GET 10287832,'NPHP1:exon_1'",df_dq_adj[(df_dq_adj['Sample']==10287832) & (df_dq_adj['gene'] == "NPHP1")])
    df_dq_adj_nowhole=df_dq_adj[df_dq_adj['adj_cnv_type'].isna()]
    df_dq_recal = df_dq_adj_nowhole[['Sample','target','dq']]
    print("remove sample+gene potentially whole gene CNV",df_dq_recal)
    df_dq_recal_back=df_dq_recal.pivot(index='Sample',columns='target')
    df_dq_recal_back=df_dq_recal_back['dq'].reset_index()
    df_dq_recal_back.columns.name=None
    df_dq_recal_back=df_dq_recal_back.set_index('Sample')
    print(df_dq_recal_back)
    df_dq_recal_back.to_csv("recal_dq.csv")
    return df_dq_recal_back

def pertargetcv_zscore_recal(df_dq_recal_back):
    print("df_dq_filterCV with dq values after removed the whole/potential gene dup/del", df_dq_recal_back)
    cols = list(df_dq_recal_back.columns)
    df_cv=pd.DataFrame()
    # for sample in df_dq_filterCV.index.values.tolist():
    for col in cols:
        cv = df_dq_recal_back[col].std(ddof=0, skipna=True) / df_dq_recal_back[col].mean(skipna=True)
        # print(sample,col,cv)
        row = {
            # 'Sample':sample,
            'target':col,
            'target_cv':cv
        }
        df_cv=df_cv.append(row,ignore_index=True)
    print(df_cv)
    df_cv_filter = df_cv[df_cv['target_cv'] <0.18]
    print(len(df_cv[df_cv['target_cv'] >=0.18])," targets failed due to cv >=0.18 !!!!",df_cv[df_cv['target_cv'] >=0.18])

    df_z = pd.DataFrame()
    df_dq_convert= df_dq_recal_back.copy()
#calculate zscore per target
    for col in cols:
        # col_zscore = col + '_zscore'
        df_z[col] = (df_dq_recal_back[col] - df_dq_recal_back[col].mean(skipna=True)) / df_dq_recal_back[col].std(ddof=0, skipna=True)
    print("df_z with zscores",df_z.head())
    df_z_convert = pd.DataFrame(index=df_z.index)
    for col in cols:
        # col_zscore = col + '_zscore'
        # col_zcnv = col +'_zcnv'
        conditions = [df_z[col] <= -3, df_z[col] >= 3,(df_z[col]>-3) & (df_z[col] <3)]
        choices = ["deletion", 'duplication','normal']
        df_z_convert[col] = np.select(conditions,choices,default=np.nan)
    print("df_z_convert with zcnv", df_z_convert.head())
    #calculate dq cv per target
    for col in cols:
        # col_dq = col+"_dqscore"
        conditions = [df_dq_recal_back[col] <= 0.65, df_dq_recal_back[col] >= 1.4, (df_dq_recal_back[col] > 0.65) & (df_dq_recal_back[col] < 1.4)]
        choices = ["deletion", 'duplication','normal']
        df_dq_convert[col] = np.select(conditions,choices,default=np.nan)
    print("df_dq_convert with dqcnv", df_dq_convert.head())
    df_dq_recal_back.reset_index(level=0, inplace=True)
    df_z.reset_index(level=0,inplace=True)
    df_z_convert.reset_index(level=0,inplace=True)
    df_dq_convert.reset_index(level=0,inplace=True)

    df_dq_filterCV2= pd.melt(df_dq_recal_back, id_vars=['Sample'], var_name='target', value_name='dq')
    df_dq_convert2=pd.melt(df_dq_convert,id_vars=['Sample'],var_name='target',value_name='dqCNV')
    df_z2=pd.melt(df_z,id_vars=['Sample'],var_name='target',value_name='zscore')
    df_z_convert2=pd.melt(df_z_convert,id_vars=['Sample'],var_name='target',value_name='zCNV')
    print(df_dq_filterCV2,df_dq_convert2,df_z2,df_z_convert2)
    df_z_res=pd.merge(df_z2,df_z_convert2,on=['Sample','target'])
    df_dq_res = pd.merge(df_dq_filterCV2,df_dq_convert2,on=['Sample','target'])
    df_z_dq = pd.merge(df_z_res,df_dq_res,on=['Sample','target'])
    df_cv_merge = df_z_dq.merge(df_cv_filter,on=['target'])

    print("Before consensus",df_cv_merge)
    df_cv_merge['consensus_CNV'] = df_cv_merge.apply(initcnv.consensus_z_dq,axis=1)
    df_cv_merge[['gene','exon']] = df_cv_merge.target.str.split(':exon_',expand=True)
    df_cv_merge['exon'] = df_cv_merge['exon'].astype(float)
    # print(set(df_z_dq['exon'].values.tolist()))
    df_cv_merge_sort=df_cv_merge.sort_values(["gene","exon"],ascending=(True,True))
    print("After consensus and sorting",df_cv_merge_sort)
    df_cv_merge_sort.to_csv("sort_recal_targetcv0.18.csv")
    return df_cv_merge_sort

def secalling(df_cv_merge_sort):
    df_consensus_nonnormal = df_cv_merge_sort[df_cv_merge_sort['consensus_CNV'] != 'normal']
    df_normal = df_cv_merge_sort[df_cv_merge_sort['consensus_CNV'] == 'normal']
    df_normal.to_csv("normal_calls.csv")
    print("cons_normal", df_normal)
    print("cons_nonnormal", df_consensus_nonnormal)
    # get potential SE list
    df_se_sub = pd.DataFrame()
    for i, row in df_consensus_nonnormal.iterrows():
        if row['target_cv'] <= 0.12 and (row['dq'] <= 0.62 or row['dq'] >= 1.45):
            df_se_sub = df_se_sub.append(row, ignore_index=True)
        elif row['target_cv'] <= 0.18 and (row['dq'] <= 0.5 or row['dq'] >= 1.5):
            df_se_sub = df_se_sub.append(row, ignore_index=True)
    print("after filter", df_se_sub)
    df_se_mark = pd.DataFrame()
    for g, df in df_se_sub.groupby(['Sample', 'gene']):
        # if g==(14076402,'FANCD2'):
        #     print(g,df)
        #     print(df_normal[(df_normal['Sample'] == 14076402) & (df_normal['gene'] == 'FANCD2') ])
        sample = g[0]
        gene = g[1]
        if len(df) == 1:  # one positive target per sample
            l, h, lb, ub = recal_SE(df_normal, sample, gene)
            df['adj_cnv_type'] = df.apply(assign_se_type, args=(l, h, lb, ub), axis=1)
            df_se_mark = df_se_mark.append(df, ignore_index=True)
        else:
            # print(g,df)
            exons_num = df['exon']
            # print(exons_num)
            exon_cl = tell_se_me(exons_num)
            # print(exon_cl)
            for i, row in df.iterrows():
                for e in exon_cl:
                    if row['exon'] in e:
                        if len(e) == 1:
                            row['adj_cnv_type'] = assign_se_type(row, l, h, lb, ub)
                        else:
                            row['adj_cnv_type'] = 'ME'
                        break
                df_se_mark = df_se_mark.append(row, ignore_index=True)
    print(df_se_mark)
    df_se_mark.to_csv("se_marked.csv")

    return df_se_mark, df_normal