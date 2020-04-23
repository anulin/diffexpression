import glob
import math
import numpy as np
import matplotlib.pyplot as plt
import os
humanfactors={}
mousefactors={}

HTOM={}
MTOH={}
n=0
with open("HOM_MouseHumanSequence.rpt.txt") as file:#соответствие генов человека и мыши
    for line in file:
        if line.split()[1]=="mouse":
            mouse=line.split()[3]
            n=line.split()[0]
        if line.split()[1] == "human" and line.split()[0]==n:
            HTOM[line.split()[3]]=mouse
            MTOH[mouse]=line.split()[3]

GoodTissues={}
n=0
humangenes={}#hgnc id в symbol
enshu={}#ensg в symbol
with open("белки.txt") as file:#hgnc id в symbol и ensg в symbol
    for line in file:
        if line.split()[1] in HTOM:
            humangenes[line.split()[0]] = line.split()[1]
            enshu[line.split()[2]]=line.split()[1]

mousegenes={}
with open("M20.txt") as file:# uniprot B symbol
    for line in file:
        if(line.split()[1]in HTOM.values()):
            mousegenes[line.split()[0]] = line.split()[1]
hcis={} #конветрация нормальных названий факторов транскрипции
mcis={}

with open("Human cistromF.txt") as file:
    for line in file:
        if(line.split()[1]in HTOM):
            hcis[line.split()[0]] = line.split()[1]
factForGen={}
for filename in hcis:
    if ("cistrome2genes_hg19\\hg19.gencode_v19." + filename + '.tsv') in glob.glob(
            os.path.join("cistrome2genes_hg19", '*.tsv')):
        with open("cistrome2genes_hg19/hg19.gencode_v19."+filename+'.tsv') as file:
            humanfactors[hcis[filename]]=set()
            for line in file:
                if(line.split()[0][:15]in enshu and int(line.split()[4])!=-1 and int(line.split()[4])<10000):
                    humanfactors[hcis[filename]].add(enshu[line.split()[0][:15]])
                    if enshu[line.split()[0][:15]] not in factForGen:
                        factForGen[enshu[line.split()[0][:15]]]={hcis[filename]}
                    else:
                        factForGen[enshu[line.split()[0][:15]]].add(hcis[filename])
with open("mouse cistromF.txt") as file:
    for line in file:
        if(line.split()[1]in HTOM.values()):
            mcis[line.split()[0]] = line.split()[1]
ensmo={}#ensmus в symbol
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
                if(line.split()[0][:18]in ensmo and int(line.split()[4])!=-1 and int(line.split()[4])<10000):
                    mousefactors[mcis[filename]].add(ensmo[line.split()[0][:18]])
                    if MTOH[ensmo[line.split()[0][:18]]] not in factForGen:
                        factForGen[MTOH[ensmo[line.split()[0][:18]]]]={MTOH[mcis[filename]]}
                    else:
                        factForGen[MTOH[ensmo[line.split()[0][:18]]]].add(MTOH[mcis[filename]])
genes={}#coдержит p-value диффэкспрессий генов
unexpressed={}
expressionH={}
expressionM={}
maxH={}
maxM={}
m=[]
h=[]
ttt={}
with open("MaxFactExprs.txt") as file:
    for line in file:
        if line.split()[0] in humanfactors:
            maxH[line.split()[0]]=float(line.split()[1])
        else:
            maxM[line.split()[0]]=float(line.split()[1])
with open("comparison.txt") as file:
    for line in file:
        if line.split()[0]=='gene':
            for i in range(len(line.split())):
                if line.split()[i][0]=='M':
                    if line.split()[i][-1]=='h':
                        h.append(i)
                    else:
                        m.append(i)
        elif line.split()[0] in humanfactors:
            unexpressed[line.split()[0]]= np.mean([float(line.split()[i]) for i in h]) < maxH[line.split()[0]] * 0.1 \
                                          and np.mean([float(line.split()[i]) for i in m]) < maxM[HTOM[line.split()[0]]] * 0.1#не экспрессируется
            expressionH[line.split()[0]]=np.mean([float(line.split()[i]) for i in h])
            expressionM[line.split()[0]]=np.mean([float(line.split()[i]) for i in m])
        if line.split()[0]!='gene' and HTOM[line.split()[0]] in mousefactors:
            ttt[line.split()[0]]=np.mean([float(line.split()[i]) for i in m]) > maxM[HTOM[line.split()[0]]] * 0.1#экспрессируется
a=0
logfC={}
with open("Card-Men_M.txt") as file:
    for line in file:
        if a!=0:
            genes[line.split()[1]]=float(line.split()[5])
            logfC[line.split()[1]]=line.split()[2]
        else:
            a=1
a=0

print(humanfactors.keys()-{MTOH[fact] for fact in mousefactors.keys()}&genes.keys())
print({MTOH[fact] for fact in mousefactors.keys()}&genes.keys()-humanfactors.keys())
st=set()#статистики
gns=0
for gn in factForGen:
    if gn in genes:
        gns+=1
    for fct in factForGen[gn]:
        if fct in genes and HTOM[fct] in mousefactors and HTOM[gn] in mousefactors[HTOM[fct]] and  ttt[fct]:
            a+=1
            st.add(fct)
print(a/gns, len(humanfactors),len(st))
print(len(mousefactors))
# sfda=[]
# for fact in st:
#     sfda.append(-math.log(genes[fact]))
# plt.hist(sfda)
# for fact in humanfactors:
#     if fact in genes:
#         sfda.append(math.log(genes[fact]))
# sfda.sort()
# plt.plot(sfda,[1-i/len(sfda)for i in range(len(sfda))])
# plt.show()
difference=[{} for i in range(10)]
st=set()
badgn=0
with open("difgenes.txt",'w') as file:
    file.write('gene\tlog(p-val)\tlogFC\tfactors/expressions;expr/max...\tAvgGeom\n')
    for gene in genes:
        if gene not in factForGen or genes[gene]>0.05:#откидываем если большое p-value гена
            continue
        file.write(gene + '\t'+str(round(math.log(genes[gene]),3))+'\t' +logfC[gene]+'\t')
        product=1
        n=0
        s=""
        for fct in factForGen[gene]:
            if (fct not in humanfactors or HTOM[fct] not in mousefactors)and fct in genes:#нет файла для фактора у человека или мышы
                st.add(fct)
                continue
            if fct in genes and( (HTOM[gene] not in mousefactors[HTOM[fct]] or gene not in humanfactors[fct])and not unexpressed[fct]):
                a = 1
                if gene in humanfactors[fct]:
                    s+=("*H*"+fct + " H:" + str(round(expressionH[fct], 3)) + ";" + str(round(expressionH[fct] / maxH[fct], 3)) + "M:" + str(
                        round(expressionM[fct], 3)) + ';' + str(round(expressionM[fct] / maxM[HTOM[fct]], 3)) + '\t'+str(
                        round(math.log(genes[fct]), 3)) + ';' + logfC[fct]+ '\t')
                else:
                    s+=('*M*'+
                        fct + " H:" + str(round(expressionH[fct], 3)) + ";" + str(round(expressionH[fct] / maxH[fct], 3)) + "M:" + str(
                            round(expressionM[fct], 3)) + ';' + str(round(expressionM[fct] / maxM[HTOM[fct]], 3)) + '\t'+str(
                        round(math.log(genes[fct]), 3)) + ';' + logfC[fct]+ '\t')
                continue #если фактор не регулирует ген у кого-то, и его експрессия велика
            elif fct in genes:
                if not unexpressed[fct]:
                    s+=(fct + " H:" + str(round(expressionH[fct], 3)) + ";" + str(
                        round(expressionH[fct] / maxH[fct], 3)) + " M:" + str(round(expressionM[fct], 3)) + ';' + str(
                        round(expressionM[fct] / maxM[HTOM[fct]], 3)) + '\t' + str(
                        round(math.log(genes[fct]), 3)) + ';' + logfC[fct]+ '\t')
                n+=1
                product*=genes[fct]
        # if a == 1 or n==0:#если есть "плохие" факторы - p-value не считаем
        #     a = 0
        #     file.write("NaN\n")
        #     continue
        product**=(1/n)
        file.write(str(product)+'\t'+s)
        for i in range(1,100,10):
            if(i/100<product):# у факторов p-value должно быть мало
                difference[(i-1)//10][gene] = fct
        file.write('\n')
print(badgn,st)
for i in range(10):
    print(difference[i].keys())


