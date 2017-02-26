import BioMod as BM
from BioMod import Seq

seqs = map(Seq.from_fasta, BM.read_fasta_file("NA.fasta"))
print [x.fasta()+"\n" for x in seqs[0:3]]
print "fasta"
octs = map(Seq.from_fasta, BM.read_fasta_file("octomers.fasta"))
entries = 0
f = open("test.fasta",'w')

def horspool(y): #This function makes the Character Shift table for each octomer
    table={}
    i=0
    x=len(y)
    for nuc in y:
        if i==x-1:
            break
        if nuc not in table.keys():
            table.update({nuc:x-1-i})
        else:
            table[nuc]=x-1-i
        i=i+1
    print table
    return table


for o in octs:
    table=horspool(o.seq)
    for s in seqs:
        i=0
        #for i in range(0,s.len() - o.len() -1,char_shift_table(o,s)):
        while i<=s.len()-o.len()-1:
            sub_seq = s.seq[i:i+o.len()]
            q=BM.get_hamming(o.seq, sub_seq, max_hamming=0)
            if  q== 'default':
                print sub_seq+"\n"+o.seq+"\n"

                #extract flanking reqions
                #f.write(str(i))
                new = Seq(o.id+"|"+s.id, s.seq[i-2:i+o.len()+2], s.chr, s.pos-5+i)
                #new = Seq(o.id+"|"+s.id, s.seq[i:i+o.len()], s.chr, s.pos)

                f.write(new.fasta())
                entries += 1
                print "Entry No: " + str(entries)
                i=i+1
            else:
                if q in table.keys():
                    i=i+table[q]
                else:
                    i=i+o.len()
