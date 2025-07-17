from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from .gtf import GTF

def create_rrna_ref(gtf: GTF):
    records=[]
    for feature in gtf.db.features_of_type('transcript'):
        seq=gtf.get_gene_seq(feature.id)
        gene_info=gtf.get_info_from_id(feature.id)
        records.append(SeqRecord(
            Seq(seq).back_transcribe(),
            id=feature.id,
            description=f"{feature.id}\t{gene_info['product'][0]}\t{gene_info['transcript_biotype'][0]}"
        ))

    with open('transcript.fna', "w") as output_handle:
        SeqIO.write(records, output_handle, "fasta")
