#author: WANG Yuanyuan
import os
import pandas as pd

def mkdir(dir):
    if not os.path.exists(dir):
        os.makedirs(dir)


def RNA_2nd(dir):
    os.chdir(dir)
    with open('rna2.txt', 'a+') as finalfile:
        list1=['chr','pos','editing_type','structure','sample']
        seq='\t'.join(list1)
        finalfile.write(seq+'\n')


        for files in os.listdir(dir):
            if os.path.isfile(files) and files.endswith('fwd.ag.rna2.fasta.txt'):
                print(files)
                name=files.split('.fwd.ag.rna2.fasta.txt')[0]
                rev=name+'.rev.tc.rna2.fasta.txt'

                with open(files,'r') as f1:
                    for i1 in f1.readlines():
                        i1=i1.strip()
                        chr=i1.split(':')[0]
                        pos=int(i1.split('\t')[0].split('-')[1])-30
                        structure=i1.split('\t')[1]
                        finalfile.write(chr+'\t'+str(pos)+'\t'+'A>G'+'\t'+structure+'\t'+name+'\n')

                with open(rev,'r') as f2:
                    for i2 in f2.readlines():
                        i2= i2.strip()
                        chr = i2.split(':')[0]
                        pos = int(i2.split('\t')[0].split('-')[1]) - 30
                        structure = i2.split('\t')[1]
                        finalfile.write(chr + '\t' + str(pos) + '\t' + 'T>C'+'\t' + structure+'\t'+name+'\n')

def reflect_to_gene(dir):
    os.chdir(dir)

    with open('list_726.txt','r') as f1:
        list1=f1.readlines()

    with open('genes_n.gff','r') as f4:
        list4=f4.readlines()

    with open('list_726_gene.txt','w') as f11:
        f11.write(list1[0].strip()+'\t'+'gene_id'+'\n')

        for i in list1[1:]:
            i=i.strip()
            chr=i.split('\t')[0]
            pos=i.split('\t')[1]
            for i4 in list4:
                i4=i4.strip()
                left=i4.split(' ')[1]
                right=i4.split(' ')[2]
                if chr == i4.split(' ')[0]:
                    if int(pos) > int(left) and int (pos) < int(right):
                        id=i4.split('geneID=')[1]

                        f11.write(i+'\t'+id+'\n')

def upsetr(dir):
    os.chdir(dir)
    df1 = pd.read_csv('nonsyn_gene_editing_level.txt', sep='\t',low_memory=False)
    df1=df1.fillna(0)
    df2 =pd.read_csv('5hosts_sample.txt',sep='\t',error_bad_lines=False)
    list1 = df2['sample'].values.tolist()
    df3=pd.DataFrame()
    df3['id']=df1['id']
    column1 = df1.columns.values.tolist()
    for i1 in list1:
        if i1 in column1:
            if i1 in df2['sample'].values.tolist():
                print (i1)
                df3[i1]=df1[i1]

    df4=df3.groupby(df3['id']).mean()
    df4.to_csv('upsetr_5hosts_average.txt', sep='\t')
    df3.to_csv('upsetr_5hosts_2.txt', sep='\t',index=None)

def change_title(dir,table):
    #基本不用改的，只要保证rnaseqlist.txt里的列名不变就可以
    os.chdir(dir)
    df1 = pd.read_csv('rnaseqlist.txt', sep='\t')
    df2 = pd.read_csv(table, sep='\t',index_col=0)
    #df2 = df2.set_index(keys=['chr'])
    df3 = pd.DataFrame(df2.T)
    df3['id'] = df3.index
    df4 = pd.merge(df1, df3, how='inner', on='id').set_index(keys=['name']).T
    df4.drop(['group','id'],axis=0, inplace=True)

    return df4

def host_specific_RES(dir):
    os.chdir(dir)

    with open('editing_level_turkey_new.txt','r') as f1:
        with open('temp.txt','w') as f2:
            dict={}
            for i in f1.readlines():
                i=i.strip()
                i=i.replace(' ','_')
                if i=='x':
                    continue
                else:
                    seq=[]
                    id=i.split(':')[-1]

                    seq.append(i.split(':')[0]+'|'+'padj|'+i.split(':')[1])

                    if id not in dict:
                        dict[id] = seq
                    else:
                        dict[id] += seq
            for k, v in dict.items():
                f2.write(k + '\t')
                seq1 = ';'.join(v)
                f2.write(seq1 + '\n')

    with open('temp.txt', 'r') as f3:
        with open('editing_level_turkey_table.txt', 'w') as f4:
            f4.write(
                'Site'+'\t'+'turkey_result' +'\t'
                'Nb_sig' + '\t' + 'Bn_sig' + '\t' + 'At_sig' + '\t' + 'Zm_sig' + '\t' + 'Ta_sig' + '\n')
            for i2 in f3.readlines():
                i2 = i2.strip()
                site=i2.split('\t')[0]
                turkey_result=i2.split('\t')[1]
                nb = str(i2.count('nb'))
                bn = str(i2.count('bn'))
                at = str(i2.count('at'))
                zm = str(i2.count('zm'))
                ta = str(i2.count('wh'))
                list2 = [nb, bn, at, zm, ta]
                seq = '\t'.join(list2)
                f4.write(site + '\t' + turkey_result + '\t' + str(seq) + '\n')

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

def merge_type_read():
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



def merge_all_needed(dir):
    #有重复值，同一位点在不同基因中出现
    os.chdir(dir)
    df1=pd.read_csv('Site_selection.txt',sep='\t')
    df2=pd.read_csv('editing_dp_ad_level.txt',sep='\t')
    df_gene=pd.read_csv('gene_editing_level.txt',sep='\t')
    df_gene=df_gene.drop_duplicates()
    df_feature=pd.read_csv('region_410.txt',sep='\t')
    df_turkey=pd.read_csv('editing_level_turkey_table.txt',sep='\t')
    df3=pd.DataFrame()
    list_col=[]
    list_col2=[]
    df_cal=pd.DataFrame()

    with open('col_name.txt','r') as f1:
        for item in f1.readlines():
            list_col.append(item.strip())

    with open('col_name2.txt', 'r') as f2:
        for item2 in f2.readlines():
            list_col2.append(item2.strip())

    list1=df2.columns.values.tolist()
    for i in list_col:
        if i in list1:
            df3[i]=df2[i]

    for i2 in list_col2:
        if i2 in list1:
            df_cal[i2]=df2[i2]

    df3.fillna('NA',inplace=True)
    #df_cal.fillna(0,inplace=True)
    # 按条件新增一列
    df3['S_NS'] = df3.missense.apply(lambda x: 'synnonymous' if x=='NA'  else 'nonsynonymous')

    df4=pd.merge(df1,df3,how='inner',on=['Site','Type'])
    df5=pd.merge(df_gene,df4,how='outer',on=['Site','Type'])
    df6=pd.merge(df_feature,df5,how='right',on=['chr','pos'])

    df7=pd.merge(df6,df_turkey,how='left',on=['Site'])
    df7.fillna('NA',inplace=True)
    df7=df7.drop_duplicates()
    df7=df7.sort_values(by='gene_id')
    df7.to_csv('merge_all_site_selection.txt',sep='\t',index=None)
    #df_cal.to_csv('pure_level.txt',sep='\t',index=None)

def cluster_res(dir):
    os.chdir(dir)
    list_result=[]
    with open('Site_selection_sort.txt','r') as f1:
        list1=f1.readlines()
    list_result.append(list1[1])
    for i1 in list1[2:]:
        i1= i1.strip()
        start=list_result[-1].split('_')[2].split('\t')[0]

        end=i1.split('_')[2].split('\t')[0].split('\t')[0]
        len1=int(end)-int(start)
        if i1.split('\t')[1] in list_result[-1] and len1< 2000:
            list_result.append('\t'+'###'+i1)
        else:
            list_result.append(i1)

    with open('site_cluster2.txt','w') as f4:
        for line in list_result[1:]:
            line=line.strip()

            if '###' in line:
                f4.write(line+'\t')
            else:
                f4.write('\n'+line)

    with open('site_cluster2.txt','r') as f5:
        list2=f5.readlines()
    with open('site_sort_cluster.txt','w') as f6:
        for i in list2[1:]:
            i=i.strip()
            count1=1+i.count('###')
            if count1>2:

                left=i.split('_')[2].split('\t')[0]
                right=i.split('_')[-1].split('\t')[0]
                chr2=i.split('_')[0]+'_'+i.split('_')[1]
                len=int(right)-int(left)
                if len<2000:
                    f6.write(chr2+'\t'+left+'\t'+right+'\t'+str(count1)+'\t'+i.split('\t')[-1]+'\n')


def main():
    """
    第一步，关于RNA二级结构的处理
    """
    dir = '/Users/yuanyuan/Downloads/yzc/RNA-editing/5hosts/'
    RNAss=dir+'2nd/'

    #RNA_2nd(RNAss)

    """
    第二步，RES该转成基因的转成基因，运行速度比较慢，20min左右，但是由于是单个位点，暂时没有更好的方法
    """
    reflect_to_gene(dir)
    """
    通用的upsetr格式转换
    """
    #upsetr(dir)
    """
    按照分组的id转换一下table
    """
    #df1='upsetr_mpmpn.txt'
    #df2='res_to_gene_editing_level.txt'

    #change_title(dir,df2).to_csv('upsetr_group2.txt', sep='\t')
    #host_specific_RES(dir)

    #merge_type_read()
    #cluster_res(dir)

    """
    需要加chro，pos，type，feature，gene_id，S_NS，AD，DP，level。
    再加上差异位点分析，如果pvalue啥的。
    type为12种，feature为CDS和5‘UTR等，S_NS为synonymous/nonsynonymous
    """
    #merge_all_needed(dir)
main()
print('finish')