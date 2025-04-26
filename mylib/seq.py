

class Aligner_FS:
    def __init__(self, pairing_rule, target_seq):
        self.pairing_rule=pairing_rule 
        self.target_seq=target_seq

    def _pairs(self, n1, n2):
        if {n1, n2} in self.pairing_rule:
            return True
        return False

    def pair(self, seq, start_i):
        results=[]
        for i in range(len(self.target_seq)):
            try:
                if self._pairs(self.target_seq[i], seq[i+start_i]):
                    results.append('+')
                else:
                    results.append('-')
            except IndexError:
                results.append('-')
        return ''.join(results)
    

if __name__=='__main__':
    from data import ASD_BA_SUB, WEAK_RNA_RULE
    aligner=Aligner_FS(WEAK_RNA_RULE, ASD_BA_SUB)
    print(aligner.pair('ZZZZZGGAGG', 0) )
