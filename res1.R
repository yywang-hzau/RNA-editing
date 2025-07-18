setwd('/Users/yuanyuan/Downloads/yzc/RNA-editing/5hosts/')


library(ggplot2)
library(stringr)
library(magrittr)
library(tidyr)


require(plyr); require(gridExtra); require(grid);

#可视化常用包，最好都加载进来
library(ggthemes)
library(Hmisc)
library(tidyverse)
library(agricolae)
library(car)
library(reshape2)
library(ggpubr)
library(dbplyr)
library(dplyr)
library(rstatix)
library(RColorBrewer)
library(ggsignif)
library(cowplot)
library(data.table)
library(ggstatsplot)

data<-read.csv('count_missense.txt',sep='\t',header=TRUE,check.names = FALSE)
print(colnames(data))
data<-set_colnames(data,c('sample','Nonsynonymous','Synonymous'))

sample<-read.csv('5hosts_sample.txt',sep='\t',header=TRUE)
print(colnames(sample))
list=unique(sample$sample)

missense<-merge(sample,data,by='sample')
missense<-data.table(missense)

missense2<-melt.data.table(missense,id.vars = c('sample','Group'), 
                           variable.name = "variants_type",value.name = 'RES')
syn<-missense2[which(missense2$variants_type=='Synonymous'),]
nonsyn<-missense2[which(missense2$variants_type=='Nonsynonymous'),]

mycolor<-brewer.pal(9,"Set1") 

#把通常用的theme封起来
theme_usual<-theme(
  panel.grid.major=element_blank(),
  panel.grid.minor=element_blank(),
  legend.title  = element_blank(),
  legend.position="right",
  plot.title = element_text(hjust = 0.5),
  text=element_text(family = 'sans',size=14),
  panel.background=element_rect(fill='transparent', color='black'),
  axis.text=element_text(color='black'),
  axis.title.x = element_blank(),
  axis.text.x=element_text(angle=0,size=14,hjust=0.5))


model<-aov(RES~Group, data=syn)
summary(model)
out <- LSD.test(model,"Group", p.adj="none")
out
#结果显示：标记字母法out$group
out$groups

#这里只需要改引号中的变量名‘Group’
sig<-out$groups
sig['Group']=row.names(sig)
sig<-sig[,-1]
mean<-out$means
mean['Group']=row.names(mean)
mean<-mean[,-1]
add<-merge(sig,mean,by='Group')
syn2<-merge(syn,add,by='Group')

model2<-aov(RES~Group, data=nonsyn)
summary(model2)
out2 <- LSD.test(model2,"Group", p.adj="none")
out2
#结果显示：标记字母法out$group
out2$groups

#这里只需要改引号中的变量名‘Group’
sig2<-out2$groups
sig2['Group']=row.names(sig2)
sig2<-sig2[,-1]
mean<-out2$means
mean['Group']=row.names(mean)
mean<-mean[,-1]
add2<-merge(sig2,mean,by='Group')
nonsyn2<-merge(nonsyn,add2,by='Group')

colourCount = length(unique(missense2$Group))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

p1<-
  ggplot(syn2,aes(x=Group,y=RES,color=Group))+
  geom_jitter(size=1) +
  geom_boxplot()+
  labs(y=expression("The number of synonymous RES"),title='')+
  theme_usual+
  scale_color_manual(values=mycolor[1:5])+theme(legend.position = '')+
  geom_text(aes(label = groups ,y= Max, x = Group),
            vjust = -0.3,size = 5,color='black')
p1

p2<-
  ggplot(nonsyn2,aes(x=Group,y=RES,color=Group))+
  geom_jitter(size=1) +
  geom_boxplot()+
  labs(y=expression("The number of nonsynonymous RES"),title='')+
  theme_usual+theme(legend.position = '')+
  scale_color_manual(values=mycolor[1:5])+
  geom_text(aes(label = groups ,y= Max, x = Group),
            vjust = -0.3,size = 5,color='black')
p2
"第一张，RES count"

plot_grid(p1,p2,labels=c('A','B'))

"region堆积柱形图"
#输入基因组区域,堆积柱形图
region<-read.csv('merge_bed_region_count.txt',sep='\t',
                 header=TRUE,check.names = FALSE)

region['intergenic']<-region$sum*region$ratio_intergenic
region<-region[,-c(7:14)]
print(colnames(region))
region<-set_colnames(region,
                     c('sample','3`UTR','5`UTR','CDS','Intron',
                       'ncRNA','Intergenic'))

region2<-data.table(merge(sample,region,by='sample')) 
region3<-melt.data.table(region2,id.vars = c('sample','Group'), 
                  variable.name = "variants_type",value.name = 'RES')%>%
  group_by(Group,variants_type) %>%
  dplyr::summarise(mean(RES)) 

#固定顺序用

r<-ggplot(region3,aes(x =Group,y = region3$`mean(RES)`,fill = variants_type))+
  geom_bar(stat='identity',position='fill')+
  scale_fill_manual(values=mycolor)+
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+
  theme_usual+labs(y='Ratio')+labs(y='Ratio')+
  theme(axis.text.x=element_text(angle=0,size=14,hjust=0.5))
r
#带百分比的饼图

region3<-set_colnames(region3,c('Group','region','mean'))
region_list<-unique(region3$region)
region_list
library(tidyverse)
library(scales)
at<-region3[which(region3$Group=='At'),]
at<-data.table(at)
#这里的rev非常重要
at$region<-factor(at$region,levels = rev(region_list))
atp<-ggplot(at, aes(x =0.5, y =mean , fill = region)) +
  geom_bar(stat = "identity",position = 'stack',width = 0.5) +
  coord_polar("y") +
  scale_fill_brewer(palette = 'Paired') +
  theme_void()+
  labs(title='At',)+
  theme(text=element_text(family = 'sans',size=14),
        legend.title=element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = '')+
  geom_text(data=at,aes(y = mean/2 +
                          c(0, cumsum(mean)[-length(mean)]),
                        label = percent(mean/ sum(mean))), size = 3)

atp


nb<-region3[which(region3$Group=='Nb'),]
#这里的rev非常重要
nb$region<-factor(nb$region,levels = rev(region_list))
nbp<-ggplot(nb, aes(x = "", y =mean , fill = region)) +
  geom_bar(stat = "identity",position = 'stack',width = 0.5) +
  coord_polar("y") + theme_void()+
  scale_fill_brewer(palette = 'Paired') +
  labs(title='Nb',)+
  theme(text=element_text(family = 'sans',size=14),
        legend.title=element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = '')+
  geom_text(data=nb,aes(y = mean/ 2 +
                          c(0, cumsum(mean)[-length(mean)]),
                        label = percent(mean/ sum(mean))), size = 3)

nbp

bn<-region3[which(region3$Group=='Bn'),]
#这里的rev非常重要
bn$region<-factor(bn$region,levels = rev(region_list))

bnp<-ggplot(bn, aes(x = "", y =mean , fill = region)) +
  geom_bar(stat = "identity",position = 'stack',width = 0.5) +
  coord_polar("y") +
  scale_fill_brewer(palette = 'Paired') +
  theme_void()+
  labs(title='Bn',)+
  theme(text=element_text(family = 'sans',size=14),
        legend.title=element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = '')+
  geom_text(data=bn,aes(y = mean/ 2 +
                          c(0, cumsum(mean)[-length(mean)]),
                        label = percent(mean/ sum(mean))), size = 3)

bnp

ta<-region3[which(region3$Group=='Ta'),]
#这里的rev非常重要
ta$region<-factor(ta$region,levels = rev(region_list))
tap<-ggplot(ta, aes(x = "", y =mean , fill = region)) +
  geom_bar(stat = "identity",position = 'stack',width = 0.5) +
  coord_polar("y") +
  scale_fill_brewer(palette = 'Paired') +
  theme_void()+
  labs(title='Ta',)+
  theme(text=element_text(family = 'sans',size=14),
        legend.title=element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = '')+
  geom_text(data=ta,aes(y = mean/ 2 +
                          c(0, cumsum(mean)[-length(mean)]),
                        label = percent(mean/ sum(mean))), size = 3)
tap

zm<-region3[which(region3$Group=='Zm'),]
#这里的rev非常重要
zm$region<-factor(zm$region,levels = rev(region_list))
zmp<-ggplot(zm, aes(x = "", y =mean , fill = region)) +
  geom_bar(stat = "identity",position = 'stack',width = 0.5) +
  coord_polar("y") +
  scale_fill_brewer(palette = 'Paired') +
  theme_void()+
  labs(title='Zm',)+
  theme(text=element_text(family = 'sans',size=14),
        legend.title=element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = '')+
  geom_text(data=zm,aes(y = mean/ 2 +
                          c(0, cumsum(mean)[-length(mean)]),
                        label = percent(mean/ sum(mean))), size = 3)

zmp

region5<-region3 %>%
  group_by(region) %>%
  dplyr::summarise(sum(mean))

region5<-set_colnames(region5,c('region','mean'))
region5$region<-factor(region5$region,levels = rev(region_list))

total<-ggplot(region5, aes(x = "", y =mean , fill = region)) +
  geom_bar(stat = "identity",position = 'stack',width = 0.5) +
  coord_polar("y") +
  scale_fill_brewer(palette = 'Paired') +
  theme_void()+
  labs(title='Total',)+
  theme(text=element_text(family = 'sans',size=14),
        legend.title=element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = '')+
  geom_text(data=region5,aes(y = mean/ 2 +
                          c(0, cumsum(mean)[-length(mean)]),
                        label = percent(mean/ sum(mean))), size = 3)
total

"第二张，关于region的饼图"
plot_grid(bnp,atp,nbp,tap,zmp,total,labels=c('A','B','C','D','E','F'))

"卡方检验"
FillNA <- function(x,n){
  x[is.na(x )]<- n;
  x
}
data<-data.table(read.csv('merge_bed_region.txt',sep='\t',
                 header=TRUE,check.names = FALSE)) %>%
  melt.data.table(id.vars = c('chr','pos','region'), 
                  variable.name = "sample",value.name = 'res') %>%
  merge(sample,by='sample') %>%
  filter(str_detect(res,'>') & !str_detect(Group,'NA'))

data$region[which(data$region=="")]<-"Intergenic"

#堆积柱形图带卡方,这里图例顺序有问题！
list_region=unique(data$region)
data$region<-factor(data$region,value=rev(list_region))

p<-ggbarstats(data, region, Group)+
  theme_usual+
  theme(legend.position = 'top')
p

'''
第三张是关于editing_type的,et是初始表格，信息最多
'''
et<-read.csv('merge_editing_type_ann.txt',sep='\t',
             header=TRUE,check.names = FALSE)
et<-data.table(et)

#这一步产生过滤后的数据，分组求和
et2<-melt.data.table(et,id.vars = c('chr','pos','missense','type'), 
                     variable.name = "sample",
                     value.name = 'RNA editing type')%>%
  filter(str_detect(`RNA editing type`,'>'))
#group_by(`RNA editing type`,type,sample) %>%
#dplyr::summarise(count=n()) %>%
#filter(count>5)

###如果要简化类型，运行下面这些行，如果不想简化，直接跑上面的melt去掉# 
et2$`RNA editing type`[which(et2$`RNA editing type`=="[F]C>T")]<-"C>T" 
et2$`RNA editing type`[which(et2$`RNA editing type`=="[R]C>T")]<-"C>T" 
et2$`RNA editing type`[which(et2$`RNA editing type`=="[F]A>T")]<-"A>T" 
et2$`RNA editing type`[which(et2$`RNA editing type`=="[R]A>T")]<-"A>T" 
et2$`RNA editing type`[which(et2$`RNA editing type`=="[F]G>T")]<-"G>T" 
et2$`RNA editing type`[which(et2$`RNA editing type`=="[R]G>T")]<-"G>T" 

et2$`RNA editing type`[which(et2$`RNA editing type`=="[F]G>A")]<-"G>A" 
et2$`RNA editing type`[which(et2$`RNA editing type`=="[R]G>A")]<-"G>A" 
et2$`RNA editing type`[which(et2$`RNA editing type`=="[F]C>A")]<-"C>A" 
et2$`RNA editing type`[which(et2$`RNA editing type`=="[R]C>A")]<-"C>A" 
et2$`RNA editing type`[which(et2$`RNA editing type`=="[F]T>A")]<-"T>A" 
et2$`RNA editing type`[which(et2$`RNA editing type`=="[R]T>A")]<-"T>A" 

et2$`RNA editing type`[which(et2$`RNA editing type`=="[F]G>C")]<-"G>C" 
et2$`RNA editing type`[which(et2$`RNA editing type`=="[R]G>C")]<-"G>C" 
et2$`RNA editing type`[which(et2$`RNA editing type`=="[F]A>C")]<-"A>C" 
et2$`RNA editing type`[which(et2$`RNA editing type`=="[R]A>C")]<-"A>C" 
et2$`RNA editing type`[which(et2$`RNA editing type`=="[F]T>C")]<-"T>C" 
et2$`RNA editing type`[which(et2$`RNA editing type`=="[R]T>C")]<-"T>C" 

et2$`RNA editing type`[which(et2$`RNA editing type`=="[F]A>G")]<-"A>G" 
et2$`RNA editing type`[which(et2$`RNA editing type`=="[R]A>G")]<-"A>G" 
et2$`RNA editing type`[which(et2$`RNA editing type`=="[F]C>G")]<-"C>G" 
et2$`RNA editing type`[which(et2$`RNA editing type`=="[R]C>G")]<-"C>G" 
et2$`RNA editing type`[which(et2$`RNA editing type`=="[F]T>G")]<-"T>G" 
et2$`RNA editing type`[which(et2$`RNA editing type`=="[R]T>G")]<-"T>G" 

et2$`RNA editing type`[which(et2$`RNA editing type`=="[F]T>A[R]A>T")]<-"T>A&A>T" 
et2$`RNA editing type`[which(et2$`RNA editing type`=="[F]G>A[R]C>T")]<-"G>A&C>T" 
et2$`RNA editing type`[which(et2$`RNA editing type`=="[F]A>G[R]T>C")]<-"A>G&T>C" 
et2$`RNA editing type`[which(et2$`RNA editing type`=="[F]T>G[R]A>C")]<-"T>G&A>C" 
et2$`RNA editing type`[which(et2$`RNA editing type`=="[F]A>T[R]T>A")]<-"T>A&A>T" 

et2$`RNA editing type`[which(et2$`RNA editing type`=="[F]A>C[R]T>G")]<-"A>C&T>G" 
et2$`RNA editing type`[which(et2$`RNA editing type`=="[F]C>A[R]G>T")]<-"C>A&G>T" 
et2$`RNA editing type`[which(et2$`RNA editing type`=="[F]C>G[R]G>C")]<-"C>G&G>C" 
et2$`RNA editing type`[which(et2$`RNA editing type`=="[F]C>T[R]G>A")]<-"C>T&G>A" 
et2$`RNA editing type`[which(et2$`RNA editing type`=="[F]G>C[R]C>G")]<-"G>C&C>G" 
et2$`RNA editing type`[which(et2$`RNA editing type`=="[F]G>T[R]C>A")]<-"G>T&C>A"
et2$`RNA editing type`[which(et2$`RNA editing type`=="[F]T>C[R]A>G")]<-"T>C&A>G"

et2<-et2 %>%
  group_by(`RNA editing type`,type,sample) %>%
  dplyr::summarise(count=n()) %>%
  filter(count>20)

et3<-merge(sample,et2,by='sample')
list_type=unique(et2$`RNA editing type`)
et3$`RNA editing type`<-factor(et3$`RNA editing type`,levels=list_type)

editing_type1<-ggboxplot(et3, x = "RNA editing type", y = "count",
                         width=0.5,
                         add = c("mean_se"),
                         color = 'Group',
                         ylab=expression('RES'),
                         xlab="")+
  geom_jitter(size=0.5,alpha=0.5)+
  facet_grid(type~.,scales = "free_y")+
  stat_compare_means(aes(group=Group),
                     label = "p.signif",
                     hide.ns=TRUE,
                     method = 'anova'
  )+
  theme_usual+
  theme(legend.position='right',
        axis.title = element_text(size=16),
        axis.text.x = element_text(angle=90,size=12,hjust=0.5)
  )+
  scale_color_brewer(palette='Set1')
editing_type1

#批量加上辅助线
len1<-length(unique(et3$`RNA editing type`))
break1<-c(seq(1,len1-1,1))
for (n in break1) 
{
  editing_type1<-editing_type1 + 
    geom_vline(xintercept = n+0.5,linetype=4,size=0.5,color='grey')
}
editing_type1

'DEGs gokegg'
见rnaediting_gokegg.r

'二级结构'

FillNA <- function(x,n){
  x[is.na(x )]<- n;
  x
}
rna2nd<-read.csv('2nd/rna2.txt',sep='\t',
                 header=TRUE,check.names = FALSE)
rna2nd$structure[which(rna2nd$structure=="h")]<-"hairpin_loop" 
rna2nd$structure[which(rna2nd$structure=="f")]<-"free" 
rna2nd$structure[which(rna2nd$structure=="i")]<-"interior_loop" 
rna2nd$structure[which(rna2nd$structure=="s")]<-"stem" 
rna2nd$structure[which(rna2nd$structure=="m")]<-"multi_loop"
rna2nd$structure[which(rna2nd$structure=="t")]<-"unknown(t)"

list_structure=unique(rna2nd$structure)

rna2nd2<-merge(sample,rna2nd,by.y='sample',all=TRUE) %>%
  group_by(structure,sample,group)  %>%
  dplyr::summarise(count=n()) %>%
  filter(!str_detect(structure,'unknown'))


#堆积柱形图带卡方

rna2nd3<-merge(sample,rna2nd,by.y='sample',all=TRUE) %>%
  filter(!str_detect(structure,'unknown'))
rna2nd3$structure<-factor(rna2nd3$structure,levels=rev(list_structure))


rnastructure2<-ggbarstats(rna2nd3, structure, Group)+
  theme_usual+labs(title='RNA secondary structure of RES')+
  scale_fill_manual(values = mycolor[c(1,2,3,4,5,9)])+
  theme(legend.position = '',
        axis.text.x=element_text(angle=0,size=12,hjust=0.5))
rnastructure2


#接下来想看一下，syn和nonsyn的二级结构差异

rna2ndsyn<-merge(et[,1:4],rna2nd3,by=c('chr','pos'),all=FALSE)
rna2ndsyn$sample<-factor(rna2ndsyn$sample,levels=list)
rna2ndsyn$structure<-factor(rna2ndsyn$structure,levels=rev(list_structure))


syn1<-ggbarstats(rna2ndsyn, structure, type, palette = 'Set1')+
  theme_usual+labs(title='RNA secondary structure of RES with synonymous or nonsynonymous mutants')+
  theme(legend.position = '',
        axis.text.x=element_text(angle=0,size=12,hjust=0.5))
syn1

rnasyn<-rna2ndsyn[which(rna2ndsyn$type=='synonymous'),]
rnanonsyn<-rna2ndsyn[which(rna2ndsyn$type=='nonsynonymous'),]
syn2<-ggbarstats(rnasyn,structure,Group, palette = 'Set1')+
  theme_usual+labs(title='RNA secondary structure of RES with synonymous mutants')+
  theme(legend.position = '',
        axis.text.x=element_text(angle=0,size=12,hjust=0.5))
syn2

syn3<-ggbarstats(rnanonsyn,structure,Group,palette = 'Set1')+
  theme_usual+labs(title='RNA secondary structure of RES with nonsynonymous mutants')+
  theme(legend.position = 'right',
        axis.text.x=element_text(angle=0,size=12,hjust=0.5))
syn3
plot_grid(rnastructure2,syn1,syn2,syn3,ncol=2,
          labels=c('A','B','C','D'))



'
PCA
'
pca<-read.csv('pca_genome.txt',sep='\t',
              header=TRUE,check.names = FALSE)
pca2<-read.csv('res_rev_pca.txt',sep='\t',
               header=TRUE,check.names = FALSE)

#pcarev<-merge(sample,pca,by='sample')
#pcafwd<-merge(sample,pca2,by='sample')

pcaplot1<-ggplot(pca,aes(x=pc1,y=pc2,color=group,fill=group))+
  geom_point(size=2)+
  geom_hline(aes(yintercept=0), colour="black", linetype="dashed")+
  geom_vline(aes(xintercept=0), colour="black", linetype="dashed")+
  theme(legend.title=element_blank(),
        legend.text = element_text(size=16),
        legend.position='right',
        legend.background = element_rect(),
        panel.background=element_rect(fill='transparent', color='gray'),
        text=element_text(family = 'sans',size=14))+
  scale_color_manual(values=mycolor)+
  labs(
    x='PC1',
    y='PC2'
  )
pcaplot1

pcaplot2<-ggplot(pca_data2,aes(x=PC1,y=PC2,color=pca2$group,fill=pca2$group))+
  geom_point(size=4)+
  geom_hline(aes(yintercept=0), colour="black", linetype="dashed")+
  geom_vline(aes(xintercept=0), colour="black", linetype="dashed")+
  theme(legend.title=element_blank(),
        legend.text = element_text(size=16),
        legend.position='right',
        legend.background = element_rect(),
        panel.background=element_rect(fill='transparent', color='gray'),
        text=element_text(family = 'sans',size=14))+
  scale_color_manual(values=mycolor)+
  labs(
    x='PC1',
    y='PC2'
  )
pcaplot2

plot_grid(pcaplot1,pcaplot2,ncol=2,labels=c('Fwd','Rev'))

