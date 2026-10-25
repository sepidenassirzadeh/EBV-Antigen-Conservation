from Bio import SeqIO

# ===== EBNA1 =====
record1 = SeqIO.read("EBV PROJECT/EBNA1.fasta", "fasta")
region1 = record1.seq[385:405]

print("EBNA1 ID:", record1.id)
print("EBNA1 region:", region1)


# ===== EBNA3A =====
record2 = SeqIO.read("EBV PROJECT/EBNA-3A.fasta", "fasta")
region2 = record2.seq[324:333]

print("EBNA3A ID:", record2.id)
print("EBNA3A region:", region2)



# ===== BZLF1 =====
record3=SeqIO.read("EBV PROJECT\BZLF1.fasta","fasta")
region3 = record3.seq[53:64]

print("BZLF1 ID:", record3.id)
print("BZLF1 region:", region3)

