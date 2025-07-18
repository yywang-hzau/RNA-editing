import os
import pandas as pd


def f(x):
    if '[F]&[R]/' in x:
        read = x.split('[F]&[R]/')[0].split('|')[0]
        level = x.split('[F]&[R]/')[0].split('|')[-1]
        dv = round(int(read) * float(level))
        return str(read) + '#' + str(dv) + '#' + str(level)
    elif '/[F]&[R]' in x:
        read = x.split('/[F]&[R]')[1].split('|')[0]
        level = x.split('/[F]&[R]')[1].split('|')[-1]
        dv = round(int(read) * float(level))
        return str(read) + '#' + str(dv) + '#' + str(level)

    elif '[F]&[R]' in x and '[F]&[R]/' not in x and '/[F]&[R]' not in x:
        left_read = x.split('[F]&[R]')[0].split('|')[0]
        right_read = x.split('[F]&[R]')[1].split('|')[0]
        left_level = x.split('[F]&[R]')[0].split('|')[-1]
        right_level = x.split('|')[-1]
        left_dv = round(int(left_read) * float(left_level))
        right_dv = round(int(right_read) * float(right_level))

        editing_level = (float(left_level) + float(right_level)) / 2
        editing_read = (int(left_read) + int(right_read)) / 2
        editing_dv = (int(left_dv) + int(right_dv)) / 2

        return str(editing_read) + '#' + str(editing_dv) + '#' + str(editing_level)
    else:
        return 'NA'+ '#' + 'NA'+ '#' + 'NA'

def merge_type():
    df0= pd.read_csv('5hosts_sample.txt', sep='\t', error_bad_lines=False)
    df1=pd.read_csv('Site_full_selection_sites_only',sep='\t')
    df2=pd.read_csv('merge_bed_ann.txt',sep='\t',low_memory=False)
    df2=df2.fillna('NA')
    df3=pd.DataFrame()
    df4=pd.DataFrame()

    list1 = df0['sample'].values.tolist()
    column1 = df2.columns.values.tolist()
    for i1 in column1[0:3]:
        df3[i1]=df2[i1]

    for i2 in list1:
        if i2 in column1[3:]:
            df_n=pd.DataFrame()
            df2[i2]=df2[i2].apply(lambda x: f(str(x)))
            df_n=df2[i2].str.split('#',expand=True)
            df3[i2+'_DP']=df_n[0]
            df3[i2+'_AD']=df_n[1]
            df3[i2+'_Level']=df_n[2]

    df3['Site']=df3['chr']+'_'+df3['pos'].map(str)
    df4 = pd.merge(df1,df3,how='inner',on='Site')

    df4.to_csv('editing_dp_ad_level.txt', sep='\t', index=None)




def main():
    os.chdir('/Users/yuanyuan/Downloads/yzc/RNA-editing/5hosts/')
    #merge_type()


main()
print('finish')