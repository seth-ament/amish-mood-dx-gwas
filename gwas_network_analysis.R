# Enrichment of GWAS risk loci for protein-protein interactions
# within a neuropsychiatry-centered gene network

# LOAD THE PACKAGES

library( biomaRt )
library( leiden )
library( igraph)

# GENE ANNOTATION

mart = useMart('ensembl')
mart = useDataset(mart=mart,'hsapiens_gene_ensembl')
anno = getBM(mart=mart,attributes=c(
  'ensembl_gene_id','ensembl_peptide_id','hgnc_symbol') )

# LOAD STRING PPI network

string = read.table('9606.protein.links.v11.5.txt.gz',header=T)
string$protein1 = gsub('9606\\.','',string$protein1)
string$protein2 = gsub('9606\\.','',string$protein2)

# LOAD gene-based p-values
# restrict to genes in the annotation

magma0 = read.table('../magma/OOA.mood_dx_broad.SCZ2.top_snp1.10kb_window.01092020.genes.out',header=T)

magma = merge( anno , magma0 , by = 1 )

string.use = string[ string[,1] %in% magma[,2] & string[,2] %in% magma[,2] , ]

magma = magma[ magma[,2] %in% union( string[,1] , string[,2] ) , ]
magma$FDR = p.adjust( magma$P )

anno.use = anno[ anno$ensembl_peptide_id %in% union(string[,1],string[,2]) , ]

# LOAD Neuropsychiatry Gene Sets

sets0 = readLines('../neuropsych_sets.txt')
sets = list()
setnames = rep(NA,length(sets0))
for( i in 1:length(sets0) ) {
  tmp = strsplit( sets0[i] , split = ' ' )[[1]]
  pept = anno.use$ensembl_peptide_id[ anno.use$ensembl_gene_id %in% tmp ]
  sets[[i]] = pept
  setnames[i] = tmp[1]
}
names(sets) = setnames

PEC = unique(unlist( sets[3:8] ))
targets = unique(unlist( sets[c(12:16,22)] ))
rare = unique(unlist(sets[c(1:2)]))
synaptome = unique(unlist(sets[9]))
sets2 = list( PEC , targets , rare , synaptome )

# which genes show up more than once

rep = table( unlist( sets2 ) )
rep2 = names(rep)[ rep > 2 ]

#### Are Amish GWAS genes enriched for interactions with known neuropsych genes?

net = rbind( string , string[,c(2,1,3)] )
amish_genes = magma$ensembl_peptide_id[ magma$P < 0.01 ]

sets[['combined']] = rep2

# Fisher's exact test

p = nGene = nEdge = or = ci.l = ci.u = rep(NA,length(sets))
for( i in 1:length(sets) ) {
  t.i = table( net[,1] %in% amish_genes , net[,2] %in% sets[[i]] )
  test = fisher.test( t.i , alternative = 'greater' )
  p[i] = test$p.value
  or[i] = test$estimate
  ci.l[i] = test$conf.int[1]
  ci.u[i] = test$conf.int[2]
  nEdge[i] = t.i[2,2]
  nGene[i] = length(sets[[i]])
  cat(i,'\n')
}
enrich = data.frame( 
  set = names(sets) , 
  nGene , 
  nEdge ,
  OR = or ,
  ci.l , ci.u , 
  P = p )
enrich = enrich[ order( enrich$P , decreasing = t ) , ]
enrich$logP = -log10(enrich$P)

write.table( enrich , quote=F , row.names=F , sep='\t' ,
	    file = 'neuropsych_list.network_enrichment.txt' )

pdf('rare.pdf')
barplot( height = enrich$logP[c(8,9,15,19)] , horiz=T )
dev.off()

pdf('gwas.pdf')
barplot( height = enrich$logP[c(2,3,10,13)] )
dev.off()

pdf('targets.pdf')
barplot( height = enrich$logP[c(5,6,20:24)] )
dev.off()

pdf('gex.pdf')
barplot( height = enrich$logP[c(7,11,12,14,16,17)] )
dev.off()

# PERMUTATION-based approach

t = table( net[,1] %in% amish_genes , net[,2] %in% rep2 )
fisher.test( t )

nperm = 100

perms0 = lapply( 1:nperm , function(i) {
   s = sample( net[,2] )
   tmp = data.frame( net[,1] , s , net[,3] )
   tmp2 = tmp[ tmp[,1] %in% amish_genes & tmp[,2] %in% rep2 , ]
   edgecount.perm = nrow(tmp2)
   tab = table( c(tmp2[,1],tmp2[,2]) )
   match = match( amish_genes , names(tab) )
   genecounts.perm = tab[ match ]
   cat(i,"\n")
   outp = list( edgecount.perm , genecounts.perm )
   return(outp)
} )

edgecount = t[2,2]

edgecounts.perm = sapply( 1:length(perms0) , 
  function(i) perms0[[i]][[1]] )
# 0 / 100 greater than observed 

fisher.test( t )$p.value
p-value = 3.483279e-19
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
 1.029893 1.046987
sample estimates:
odds ratio 
  1.038413 

#### Which Amish GWAS genes are central in the neuropsych net?
# centrality analysis
# export top 250 genes as network for Cytoscape

psychnet2 = string[ string[,1] %in% rep2 | string[,2] %in% rep2 , ]
psychnet2 <- graph_from_data_frame(psychnet2 , directed = FALSE)
 
cent2 = eigen_centrality( psychnet2 )$vector
cent2 = data.frame( peptide = names(cent2) ,
                  centrality = cent2 )
magma2 = merge( cent2 , magma , by.x = 1 , by.y = 2 , all.y = 2 )
magma2$centrank = rank(1-magma2$centrality)
magma2$prank = rank(magma2$P)
magma2$medrank = apply( magma2[,c('centrank','prank')] , 1 , median )
magma2$rank = rank( magma2$medrank )
magma2 = magma2[ order( magma2$medrank , decreasing = F ) , ]
 
magma2.sig = magma2[ magma2$P < 0.01 , ]
magma2.sig = magma2.sig[ duplicated( magma2.sig$hgnc_symbol ) == F &
			magma2.sig$hgnc_symbol != '' , ]
# magma2.sig[ 1:250 ,c(4,5,6,2,12)]
topgenes = magma2.sig$peptide[ 1:250 ]  
 
magma3 = merge( magma2[,c(3,1,2,4,14:17)] , magma0 , by = 1 , all.y = T )
magma3 = magma3[ order( magma3$P ) , ]
magma3$top250 = magma3$peptide %in% topgenes
magma3 = magma3[ duplicated( magma3$ensembl_gene_id ) == F , ]

write.table( magma3 , row.names=F , quote=F , sep='\t' ,
	     file = 'OOA.mood_dx_broad.SCZ2.top_snp1.10kb_window.01092020.with_gene_anno_and_centrality.genes.out.txt' )




g = string[ string[,1] %in% topgenes & string[,2] %in% topgenes ,] 
g = graph_from_data_frame( g , directed = F )

leid = leiden( g , resolution_parameter = 1 )

df.leid = data.frame( peptide = names(V(g)) , cluster = leid )
df.leid = merge( df.leid , magma2.sig , by = 'peptide' )

write.table( df.leid , quote=F , row.names=F , sep='\t' ,
	    file = 'OOA.mood.string.top250.nodes.txt' )

g = string[ string[,1] %in% topgenes & string[,2] %in% topgenes ,] 
g = merge( g , anno.use , by.x = 1 , by.y = 2 )
g = merge( g , anno.use , by.x = 2 , by.y = 2 )
g = g[ , c(5,7,3) ]
colnames(g) = c('node1','node2','combined_score')

write.table( g , quote=F , sep='\t' , row.names=F , 
	    file = 'OOA.mood.string.top250.edges.txt' )





