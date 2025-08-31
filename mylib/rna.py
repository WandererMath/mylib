import RNA
import ray


# Visualization tool: http://rna.tbi.univie.ac.at/forna/

@ray.remote
def seq_to_structure_energy(seq):
    structure, mfe = RNA.fold(seq)  # predict MFE structure for one sequence
    return structure, mfe

if __name__ == "__main__":
    seq1 = "CCCCAAAGGGGAAAAAAAAAAAAAAAAGGG"
    seq2 = "AAAACCCUUUUUUUUUUUUUU"
    seq = seq1 + "&" + seq2

    # Use cofold directly
    structure, mfe = RNA.cofold(seq)
    structure = structure[:len(seq1)] + "&" + structure[len(seq1):]
    print("Sequence:", seq)
    print("Structure:", structure)
    print("MFE (kcal/mol):", mfe)
