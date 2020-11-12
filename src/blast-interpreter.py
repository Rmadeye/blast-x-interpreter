from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

protein_query = "MAMLQTNLGFITSPTFLCPKLKVKLNSYLWFSYRSQVQKLDFSKRVNRSYKRDALLLSIKCSSSTGFDNSNVVVKEKSVSVILLAGGQGKRMKMSMPKQYIPLLGQPIALYSFFTFSRMPEVKEIVVVCDPFFRDIFEEYEESIDVDLRFAIPGKERQDSVYSGLQEIDVNSELVCIHDSARPLVNTEDVEKVLKDGSAVGAAVLGVPAKATIKEVNSDSLVVKTLDRKTLWEMQTPQVIKPELLKKGFELVKSEGLEVTDDVSIVEYLKHPVYVSQGSYTNIKVTTPDDLLLAERILSEDS"

nucleotide_query = "ATGGCGATGCTTCAGACGAATCTTGGCTTCATTACTTCTCCGACATTTCTGTGTCCGAAGCTTAAAGTCAAATTGAACTCTTATCTGTGGTTTAGCTATCGTTCTCAAGTTCAAAAACTGGATTTTTCGAAAAGGGTTAATAGAAGCTACAAAAGAGATGCTTTATTATTGTCAATCAAGTGTTCTTCATCGACTGGATTTGATAATAGCAATGTTGTTGTGAAGGAGAAGAGTGTATCTGTGATTCTTTTAGCTGGAGGTCAAGGCAAGAGAATGAAAATGAGTATGCCAAAGCAGTACATACCACTTCTTGGTCAGCCAATTGCTTTGTATAGCTTTTTCACGTTTTCACGTATGCCTGAAGTGAAGGAAATTGTAGTTGTATGTGATCCTTTTTTCAGAGACATTTTTGAAGAATACGAAGAATCAATTGATGTTGATCTTAGATTCGCTATTCCTGGCAAAGAAAGACAAGATTCTGTTTACAGTGGACTTCAGGAAATCGATGTGAACTCTGAGCTTGTTTGTATCCACGACTCTGCCCGACCATTGGTGAATACTGAAGATGTCGAGAAGGTCCTTAAAGATGGTTCCGCGGTTGGAGCAGCTGTACTTGGTGTTCCTGCTAAAGCTACAATCAAAGAGGTCAATTCTGATTCGCTTGTGGTGAAAACTCTCGACAGAAAAACCCTATGGGAAATGCAGACACCACAGGTGATCAAACCAGAGCTATTGAAAAAGGGTTTCGAGCTTGTAAAAAGTGAAGGTCTAGAGGTAACAGATGACGTTTCGATTGTTGAATACCTCAAGCATCCAGTTTATGTCTCTCAAGGATCTTATACAAACATCAAGGTTACAACACCTGATGATTTACTGCTTGCTGAGAGAATCTTGAGCGAGGACTCATGA"

def blast_input():
    my_query = Seq(protein_query)
    # query = SeqIO.read("sample.fasta", format="fasta")
    try:
        print("Query sent to NCBIWWW against swissprot database")
        # print(type(my_query))
        # print(type(query.seq))
        result_handle = NCBIWWW.qblast("blastp", "swissprot", my_query)
        blast_result = open("my_blast_result.xml", "w")
        blast_result.write(result_handle.read())
        blast_result.close()
        result_handle.close()
    except ValueError:
        my_translated_query = my_query.translate()
        print(f"DNA/RNA sequence detected. Sequence translated to : {my_translated_query}")
        result_handle = NCBIWWW.qblast("blastp", "swissprot", my_translated_query)
        blast_result = open("my_blast_result.xml", "w")
        blast_result.write(result_handle.read())
        blast_result.close()
        result_handle.close()


def read_blast_output():
    result = open("my_blast_result.xml")
    blast_records = NCBIXML.parse(result)
    number_of_alignments = 0
    hit_ids = []
    for rec in blast_records:
        for alignment in rec.alignments:
            hit_ids.append(alignment.accession)
            number_of_alignments += 1
            if number_of_alignments == 3: # Take three top alignments and stop
                break
    print(f"Top hits: {hit_ids}")
    return hit_ids

def read_target_db(input_list: list):
    with open('../target', 'r') as df:
        for line in df:
            line = line.split()
            try:
                for entry in input_list:
                    if entry in line[2]:
                        return print(f"Entry {entry} present in database!")
            except:
                pass

if __name__ == '__main__':
    blast_input()
    read_target_db(read_blast_output())
