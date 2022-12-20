from Bio import pairwise2
from Bio.pairwise2 import format_alignment 
from Bio.Seq import Seq

def alinhamento(seqt, numseq):
    alpha = 1
    beta = 0
    delta = -1

    score = 0
    scoreaux = 0
    count1 = 0
    count2 = 0
    scoretotal = 0

    n = numseq     # numero de sequencias
    k = 0
    scoretotal = 0

    while n > 0:
        if n == numseq:  # aumenta-se 0.5 quando gaps estao junto para cumprir o requisito
            alignments = pairwise2.align.globalms(seqt[k], seqt[k + 1], alpha, beta, delta,0)
            for alignment in alignments:
                #substituindo pelas sequencias alinhadas
                if count1 == 0:
                    seqt[k] = alignment.seqA
                    seqt[k + 1] = alignment.seqB
                    score = alignment.score
                    count1 = count1 + 1
                else:
                    if alignment.score >= score:
                        seqt[k] = alignment.seqA
                        seqt[k + 1] = alignment.seqB
                        score = alignment.score
            n = n - 2
            k = k + 1
            #calculo score total
            scoretotal = scoretotal + alignment.score
        else:
            alignments = pairwise2.align.globalms(seqt[k], seqt[k + 1], alpha, beta, delta,0)        
            for alignment in alignments: 
                if count2 == 0:
                    seqt[k] = alignment.seqA
                    seqt[k + 1] = alignment.seqB
                    scoreaux = alignment.score
                    count2 = count2 + 1
                else:    
                    if alignment.score >= scoreaux:   
                        seqt[k] = alignment.seqA
                        seqt[k + 1] = alignment.seqB
                        scoreaux = alignment.score    
            n = n - 2
            k = k + 1
            scoretotal = scoretotal + alignment.score 

    #guardar o valor total do score
    scoretotal = "Score Total: " + str(scoretotal)        
    seqt.append(scoretotal)
    return seqt