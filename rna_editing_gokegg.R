'''
从RES中获得各类信息中获取相应的基因列表，进行功能富集GO,KEGG。
同时根据数据分组尝试各类可视化
'''
"""
93行以前是配置GO,KEGG的文件，文件夹中的go_v1.0.txt文件是之前的go注释文件，富集效果不理想
2.0是新比对了20年NR库之后，根据idmapping_selected.tab.gz提取的
提取使用的是UniProt2GO_annotate.py，脚本在文件夹中
"""
#切换工作路径
setwd('/Users/yuanyuan/Downloads/yzc/RNA-editing/5hosts/')

#go注释
library(AnnotationForge)
library(clusterProfiler)
library(enrichplot)
#下面三个是做kegg接口的
library(stringr)
library(magrittr)
library(tidyverse)

library(data.table)
library(tidyr)
#venn图
library(UpSetR)
require(plyr); require(gridExtra); require(grid);
library(Vennerable)

#可视化常用包，最好都加载进来
library(ggplot2)
library(ggthemes)
library(Hmisc)
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

#做go接口

fChr <- read.csv('gokegg/fchr.txt',sep='\t',header = TRUE)
fGO <- read.csv('gokegg/fgo_2.0.txt',sep='\t',header=TRUE)

## Now make a list
## data <- list(gene_info=fSym, chromosome=fChr, go=fGO)
genus <- "Myzus"
species <- "per"
dbName <- AnnotationForge:::.generateOrgDbName(genus,species)
## this becomes the file name for the DB
dbfile <- paste(dbName, ".sqlite", sep="")

makeOrgPackage( chromosome=fChr, go=fGO,
                version="0.1",
                maintainer="yy<so@someplace.org>",
                author="yy<so@someplace.org>",
                outputDir = ".",
                tax_id="13164",
                genus=genus,
                species=species,
                goTable="go")

## then you can install on the return value
install.packages("./org.Mper.eg.db",repos = NULL,type = 'source')

"""
如果本地已经有org.Mper.eg.db，
可以直接library org.Mper.eg.db不需要再次注释！
"""
library(org.Mper.eg.db)
keytypes(org.Mper.eg.db)


#kegg通路相对容易些，但是每次都要倒入
gene_kegg <- read.csv('/Users/yuanyuan/Downloads/cyz/RNA-editing/mpmpn2/gokegg/mper-kegg.txt',
                      sep=" ",check.names=F) %>%
  dplyr::select(gene,name) %>%
  filter(str_detect(name,"K")) %>%
  set_colnames(c('GENE_ID','kegg'))
#转换k number和pathway
kegg_path<-bitr_kegg(gene_kegg$kegg, "kegg", "Path", "ko")
#合并一下并按照path.gene-id的默认格式排列，列的顺序必须是这样子
gene_pathway<-merge(gene_kegg,kegg_path,by='kegg')%>%
  dplyr::select(3,2)

pathway_name <- clusterProfiler:::kegg_list("pathway") %>%
  mutate(across("from",str_replace,"path:map","ko")) %>% 
  set_colnames(c("path_id","path_name"))

'''
接下来做一些尝试，第一个是单独一列基因列表,不管分组
'''

re1<-read.csv('specific_genelist.txt',sep='\t',check.names=FALSE)
#先录入数据，再设置行名，这样可以保持第一列还是原始值
row.names(re1)<-re1[,1]


#以这一列前2000个基因为例子
x<-re1$MYZPE13164_O_EIv2.1_0002220

#挑选出的1000个基因的go分组情况
ggo <- groupGO(x, OrgDb = 'org.Mper.eg.db',
               ont='BP',
              level=3,
               keyType = 'GID')
head(ggo)

#GO富集分析
ego <- enrichGO(x, 
                OrgDb = 'org.Mper.eg.db', 
                ont='ALL',pAdjustMethod = 'BH',
                pvalueCutoff = 0.05, 
               qvalueCutoff = 0.2,
               keyType = 'GID')
#看一下是否有结果
head(ego)


#dotplot,4.0以后版本是需要，library(enrichplot),这里的10可以调整，但是太多会很密集。
#分面
dotplot(ego,title = '', 
        showCategory = 10,
        color = 'p.adjust', 
        split='ONTOLOGY')+ 
  facet_grid(ONTOLOGY~.,scale="free")

#不分面
dotplot(ego,title = '', 
        showCategory = 10,
        color = 'p.adjust', 
        split='ONTOLOGY')

#柱形图，以count为x轴
barplot(ego, title = '', 
        showCategory = 10)

#柱形图，以q value 为x轴
mutate(ego, qscore = -log(p.adjust, base=10)) %>% 
  barplot(x="qscore")

#如果想知道go结果中基因的关联
library(ggupset)
upsetplot(ego)
#如果想保留go结果也可以导出
write.table(ego@result,file = "ego.txt",sep = "\t",
            quote = FALSE, append = FALSE, na = "NA")

#kegg富集
kegg <- enricher(x,
                 TERM2GENE = gene_pathway,
                 TERM2NAME = pathway_name,
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.05)

#kegg点图
dotplot(kegg, showCategory=30)

#kegg柱形图，以count为x轴
barplot(kegg, title = '', 
        showCategory = 10)

#如果想导出kegg的结果
write.table(kegg@result,file = "commondegs-kegg.txt",
            sep = "\t",quote = FALSE, append = FALSE, na = "NA")



'''
开始准备分组情况的,go/kegg分析
首先第一步是将示例数据中，每一列一个寄主的go/kegg写入文件中
然后对带有分组标签的go/kegg@result，再进行做图,思路和常规可视化比较像了
'''

data<-read.csv('vsbn_degs.txt',
               sep="\t",check.names=F)


for (n in c(1:length(data)))
{
  genelist<-data[,n]
  go <- enrichGO(genelist, OrgDb = 'org.Mper.eg.db', 
                 ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05, 
                 qvalueCutoff = 0.2,keyType = 'GID')
  gore<-as.data.frame(go@result,row.names=go@result$ID)
  if (nrow(gore)!=0) 
    {
    gore['label']<-colnames(data)[n]
    write.table(gore,file = "go_5host.txt",sep = "\t",quote = FALSE,
                append = TRUE, na = "NA",col.names = FALSE)
  }
  
  kegg <- enricher(genelist,
                   TERM2GENE = gene_pathway,
                   TERM2NAME = pathway_name,
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05)
  kere<-as.data.frame(kegg@result,row.names = kegg@result$ID)
  if (nrow(kere)!=0) 
    {
    kere['label']<-colnames(data)[n]
    write.table(kere,file = "kegg_5host.txt",sep = "\t",
                quote = FALSE, append = TRUE, na = "NA",col.names = FALSE)
  }
}


###kegg作图,由于kegg这里冗余的信息比较多，所有库太大，导致qvalue过大，因此暂时用pvalue
colkegg<-c('id','ID','Description','GeneRatio',
       'BgRatio',	'pvalue',	'padjust','qvalue',
       'geneID','Count','label')

hostskegg<-read.csv('kegg_5host.txt',sep='\t',
                    header = FALSE,col.names = colkegg)

hostskegg2<-hostskegg[which(hostskegg$qvalue <0.05), ]

keggplot<-ggplot(hostskegg2, aes(x=label, y=Description)) +
  geom_point(aes(size=Count,color=qvalue)) +
  scale_colour_gradient(low='blue',high='red')+
  theme_bw()+labs(
    size="Count",
    x="",
    y="",
    title="KEGG terms for DEGs among 5 hosts")+
  theme(
    plot.title = element_text(hjust = 0.5),
    text=element_text(family = 'sans',size=14),
    legend.position='right',
    panel.background=element_rect(fill='transparent', color='gray'),
    axis.text=element_text(color='black')
  )
keggplot

###go作图,这里多加一个分面吧，不需要可以去掉,换一种红绿色试试看
colgo<-c('id','ONTOLOGY','ID','Description','GeneRatio',
           'BgRatio',	'pvalue',	'padjust','qvalue',
           'geneID','Count','label')

hostsgo<-read.csv('go_5host.txt',sep='\t',
                    header = FALSE,col.names = colgo)

hostsgo2<-hostsgo[which(hostsgo$pvalue <0.05), ]

goplot<-ggplot(hostsgo2, aes(x=label, y=Description)) +
  geom_point(aes(size=Count,color=qvalue)) +
  scale_colour_gradient(low='green',high='red')+
  facet_grid(ONTOLOGY~.,scale="free")+
  theme_bw()+labs(
    size="Count",
    x="",
    y="",
    title="GO terms for SNPs of 9 hosts")+
  theme(
    plot.title = element_text(hjust = 0.5),
    text=element_text(family = 'sans',size=14),
    legend.position='right',
    panel.background=element_rect(fill='transparent', color='gray'),
    axis.text=element_text(color='black')
  )
goplot

'''
其实y叔的新包里有“Biological theme comparison”这一部分，实现起来比上面容易很多，
再来测试一下吧
其实结果和上面的方法有大部分交集，但是compareCluster很多参数无法设定，比如
ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2
结果和上面的方法不完全相同
'''
#录入细表，格式稍微不同
data2<-read.csv('9hosts_2.txt',sep='\t',header=TRUE)
#注意这里只有go,kegg不行的，因为非模式物种只能做enricher,enrichKEGG不行的哦
formula_res <- compareCluster(Id~species, data=data2, 
                              fun="enrichGO",
                              OrgDb = 'org.Mper.eg.db',
                              keyType = 'GID')
head(formula_res)

#因为这里只有寄主一个分组，如果有其他的分组也可以分面哦，比如上调下调
dotplot(formula_res) 

'''
upsetr
'''
up<-read.csv('upsetr_group.txt',sep='\t',header=TRUE,row.names = 1)

#做一个简单的venn图
upset(up, nsets=30,nintersects = 70,
      keep.order=F,
      mb.ratio = c(0.5, 0.5),
      line.size = 1,
      order.by = c('degree','freq'), decreasing = c(TRUE,TRUE),
      text.scale = 1,show.numbers = 'yes')
#这部分交集没有特别值得做的
'''
试一下根据mpmpn分组的数据,没有富集结果
'''
mpmpn<-read.csv('upsetr_group2.txt',sep='\t',header=TRUE,check.names = FALSE)
list<-read.csv('rnaseqlist.txt',sep='\t',header=TRUE)
mpmpn<-data.table(mpmpn) 
mpmpn2<- mpmpn%>%
  melt.data.table(mpmpn,id.vars = 1, 
                  measure.vars = 2:13,
                  variable.name = "name",
                  value.name = 'RNA editing level') %>%
  filter(`RNA editing level` >0.2)
  mpmpn3<-merge(mpmpn2,list,by='name')
  mpmpn3<-mpmpn3[order(mpmpn3$`RNA editing level`,decreasing = T),]
  mpmpn3<-mpmpn3[!duplicated(mpmpn3[,c('V1','group')]),]
  
mp<-mpmpn3$V1[which(mpmpn3$group=="Host_Mp_2021")]
mpn<-mpmpn3$V1[which(mpmpn3$group=="Host_Mpn_2021")]
venn1<-Venn(list('Mp'=mp,'Mpn'=mpn))
plot(venn1,doWeight=T)
in1<-venn1@IntersectionSets[["11"]]
in2<-venn1@IntersectionSets[["10"]]
in3<-venn1@IntersectionSets[["01"]]
in3go <- enrichGO(in3, 
                OrgDb = 'org.Mper.eg.db', 
                ont='ALL',pAdjustMethod = 'BH',
                pvalueCutoff = 0.05, 
                qvalueCutoff = 0.2,
                keyType = 'GID')
head(in3go)

kegg <- enricher(in2,
                 TERM2GENE = gene_pathway,
                 TERM2NAME = pathway_name,
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.05)
head(kegg)




  res2 <- compareCluster(V1~group, data=mpmpn3, 
                                fun="enrichGO",
                                OrgDb = 'org.Mper.eg.db',
                                keyType = 'GID')
  
  head(res2)
  dotplot(res2)  
  
'''
试一下gesGO,流程好了，但数据集不合适
'''
#这里注意一下，上面是data.table过的，所以没有办法添加行名
mpmpn3<-data.frame(mpmpn3)
row.names(mpmpn3)<-mpmpn3$V1

mpmpn4<-as.numeric(mpmpn3$RNA.editing.level)
names(mpmpn4)<-mpmpn3$V1
class(mpmpn4)

ego3 <- gseGO(geneList     = mpmpn4,
              
              OrgDb        = org.Mper.eg.db,
              ont          = "CC",
              minGSSize    = 50,
              maxGSSize    = 500,
              pvalueCutoff = 0.5,
              keyType = 'GID',
              verbose      = FALSE)
