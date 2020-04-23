import glob
import os
import numpy as np

humanfactors={}
mousefactors={}
mtissues={}
htissues={}
factors=set()
HTOM={}
n=0

with open("HOM_MouseHumanSequence.rpt.txt") as file:#соответствие генов человека и мыши
    for line in file:
        if line.split()[1]=="mouse":
            mouse=line.split()[3]
            n=line.split()[0]
        if line.split()[1] == "human" and line.split()[0]==n:
            HTOM[line.split()[3]]=mouse
with open("факторы с id.txt") as file:
    for line in file:
        if(line.split()[0] in HTOM):
            factors.add(line.split()[0])
            factors.add(HTOM[line.split()[0]])
GoodTissues={}
n=0
name=""
with open("куски.txt") as file:
    for line in file:
        if n==1:
            name=line.split()[0]
            for t in line.split()[1:]:
                mtissues[t]=line.split()[0]
            n=0
        else:
            n=1
            for t in line.split()[0:]:
                htissues[t]=name
humangenes={}
enshu={}
with open("белки.txt") as file:#hgnc id в symbol и symbol в ensg
    for line in file:
        if line.split()[1] in HTOM:
            humangenes[line.split()[0]] = line.split()[1]
            enshu[line.split()[2]]=line.split()[1]

mousegenes={}
with open("M20.txt") as file:# uniprot B symbol
    for line in file:
        if(line.split()[1]in HTOM.values()):
            mousegenes[line.split()[0]] = line.split()[1]
hcis={}#конветрация нормальных названий факторов транскрипции
mcis={}
with open("Human cistromF.txt") as file:
    for line in file:
        if(line.split()[1]in HTOM):
            hcis[line.split()[0]] = line.split()[1]
for filename in hcis:
    if ("cistrome2genes_hg19\\hg19.gencode_v19." + filename + '.tsv') in glob.glob(
            os.path.join("cistrome2genes_hg19", '*.tsv')):
        with open("cistrome2genes_hg19/hg19.gencode_v19."+filename+'.tsv') as file:
            humanfactors[hcis[filename]]=set()
            for line in file:
                if(line.split()[0][:15]in enshu):
                    humanfactors[hcis[filename]].add(enshu[line.split()[0][:15]])
with open("mouse cistromF.txt") as file:
    for line in file:
        if(line.split()[1]in HTOM.values()):
            mcis[line.split()[0]] = line.split()[1]
ensmo={}
with open("белкиMous.txt") as file:
    for line in file:
        if ("uniprot:"+line.split()[0]) in mousegenes:
            ensmo[line.split()[1]]=mousegenes["uniprot:"+line.split()[0]]
for filename in mcis :
    if ("cistrome2genes_mm9\\mm9.gencode_vM1." + filename + '.tsv') in glob.glob(
            os.path.join("cistrome2genes_mm9", '*.tsv')):
        with open("cistrome2genes_mm9/mm9.gencode_vM1."+filename+'.tsv') as file:
            mousefactors[mcis[filename]]=set()
            for line in file:
                if(line.split()[0][:18]in ensmo):
                    mousefactors[mcis[filename]].add(ensmo[line.split()[0][:18]])
exprf={}
expr={}
positions=[]
names=[]
hgcids=set()
ash=set()
with open("gene_expressions.tsv") as file:
    i=0
    for line in file:
        while i< len(line.split()):
            if len(line.split()[i].split("."))>3 and line.split()[i].split(".")[3]in htissues:# если есть в тканяx
                positions.append(i)
                names.append(htissues[line.split()[i].split(".")[3]])
            i+=1
        if line.split()[3] in humangenes and humangenes[line.split()[3]] in humanfactors and(# при повторе максимизируется средняя экспрессия только по исследуемым тканям!
                humangenes[line.split()[3]] not in exprf or np.mean(exprf[humangenes[line.split()[3]]])<np.mean([float(line.split()[j]) for j in positions])):
            if humangenes[line.split()[3]] in exprf:
                ash.add(humangenes[line.split()[3]])
            exprf[humangenes[line.split()[3]]] = []
            for j in positions:
                exprf[humangenes[line.split()[3]]].append(float(line.split()[j]))
        elif line.split()[3] in humangenes and humangenes[line.split()[3]] in enshu.values()and (
                humangenes[line.split()[3]] not in expr or np.mean(expr[humangenes[line.split()[3]]])<np.mean([float(line.split()[j]) for j in positions])):
            if humangenes[line.split()[3]] in expr:
                ash.add(humangenes[line.split()[3]])
            expr[humangenes[line.split()[3]]] = []
            for j in positions:
                expr[humangenes[line.split()[3]]].append(float(line.split()[j]))
positionsm=[]
namesm=[]
exprmf={}
exprm={}
asd=set()
with open("gene_expressionsM.tsv") as file:
    i = 0
    for line in file:
        while i < len(line.split()):
            if len(line.split()[i].split("."))>3 and line.split()[i].split(".")[3] in mtissues:
                positionsm.append(i)
                namesm.append(mtissues[line.split()[i].split(".")[3]])
            i += 1
        if line.split()[3] in mousegenes and mousegenes[line.split()[3]] in mousefactors and (# при повторе максимизируется средняя экспрессия только по исследуемым тканям!
                mousegenes[line.split()[3]] not in exprmf or np.mean(exprmf[mousegenes[line.split()[3]]])<np.mean([float(line.split()[j]) for j in positionsm])):
            if mousegenes[line.split()[3]] in exprmf:
                asd.add(mousegenes[line.split()[3]] )
            exprmf[mousegenes[line.split()[3]]] = []
            for j in positionsm:
                exprmf[mousegenes[line.split()[3]]].append(float(line.split()[j]))
        elif line.split()[3] in mousegenes and mousegenes[line.split()[3]] in ensmo.values() and(
                mousegenes[line.split()[3]] not in exprm or np.mean(exprm[mousegenes[line.split()[3]]])<np.mean([float(line.split()[j]) for j in positionsm])):
            if mousegenes[line.split()[3]] in exprm:
                asd.add(mousegenes[line.split()[3]] )
            exprm[mousegenes[line.split()[3]]]=[]
            for j in positionsm:
                exprm[mousegenes[line.split()[3]]].append(float(line.split()[j]))

mth={}
i=0
print(len(asd),len(ash))
# while i>-1:
#     i+=1
while i<len(names):
    j=0
    while j<len(namesm):
        if names[i]==namesm[j]:
            mth[i]=j
        j+=1
    i+=1
results=[]
fcounter=set()
badfcounter=set()
genes=set()
for fact in exprf:
    if HTOM[fact]in exprmf or HTOM[fact] in exprm:
        for i in range(1):# для каждой ткани
            #results.append(["#", names[i], fact, exprf[fact][i], exprmf[HTOM[fact]][mth[i]]])
            if """float(exprmf[HTOM[fact]][mth[i]]) == float(exprf[fact][i]) or float(
                    exprmf[HTOM[fact]][mth[i]]) != 0 and float(exprf[fact][i]) / float(
                    exprmf[HTOM[fact]][mth[i]]) < 1.2 and float(exprf[fact][i]) / float(
                    exprmf[HTOM[fact]][mth[i]]) > 0.8""" :
                if i==0:
                    fcounter.add(fact)
                for gene in humanfactors[fact]:
                    if gene in expr and HTOM[gene] in exprm and( """float(exprm[HTOM[gene]][mth[i]]) == float(expr[gene][i]) or float(
                            exprm[HTOM[gene]][mth[i]]) != 0 and (float(expr[gene][i]) / float(
                            exprm[HTOM[gene]][mth[i]]) >= 5 or float(expr[gene][i]) / float(
                            exprm[HTOM[gene]][mth[i]]) <= 0.2)"""):
                        genes.add(gene)
                        "results.append([fact, names[i], gene, expr[gene][i], exprm[HTOM[gene]][mth[i]]])"
                    elif (gene in expr or gene in exprf) and (HTOM[gene] in exprm or HTOM[gene] in exprmf):
                        fcounter.add(gene)
with open("comparison.txt",'w') as file:
    file.write('gene ')
    for i in names:
        file.write(i+'_h ')
    for i in namesm:
        file.write(i+'_m ')
    file.write('\n')
    for gene in genes:
        if gene in fcounter:
            print(';;;;;;;;;;;;;;;;;;;;;;;',gene)
            continue
        file.write(gene+" ")
        for i in range(len(names)):
            file.write(str(expr[gene][i]) +' ')
        for i in range(len(namesm)):
            file.write(str(exprm[HTOM[gene]][i])+' ')
        file.write('\n')
    for gene in fcounter:
        mark=''
        file.write(gene + " ")
        if gene in exprf:
            for i in range(len(names)):
                file.write(str(exprf[gene][i]) + ' ')
        else:
            mark+='M_'
            for i in range(len(names)):
                file.write(str(expr[gene][i]) + ' ')
        if HTOM[gene] in exprmf:
            for i in range(len(namesm)):
                file.write(str(exprmf[HTOM[gene]][i]) + ' ')
        else:
            mark+='H_'
            for i in range(len(namesm)):
                file.write(str(exprm[HTOM[gene]][i]) + ' ')
        file.write(mark+'\n')
print(len(fcounter), fcounter)