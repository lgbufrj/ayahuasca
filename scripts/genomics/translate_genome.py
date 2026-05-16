from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import seq3

# Função para traduzir a sequência de DNA para proteína
def translate_sequence(dna_seq):
    dna_seq = Seq(dna_seq)
    return dna_seq.translate()

translated_records = []

# Lê a sequência de DNA do arquivo FASTA, traduz e armazena as sequências traduzidas
for seq_record in SeqIO.parse("cds.fasta", "fasta"):
    ID = seq_record.description
    dna_seq = str(seq_record.seq)

    aa_length = len(dna_seq) % 3
    if aa_length != 0:
        dna_seq += 'N' * (3 - aa_length)
        
    aa_1letter = translate_sequence(dna_seq)

    #Cria um novo SeqRecord para a sequência de proteína traduzida
    new_record = SeqIO.SeqRecord(Seq(str(aa_1letter)), id=seq_record.id, description=ID)
        
    translated_records.append(new_record)

SeqIO.write(translated_records, "prot.fasta", "fasta")
