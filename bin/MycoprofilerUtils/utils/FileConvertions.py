
def Fasta_to_Phylip(fasta, outFile):
    '''
        Converts Fasta to Phylip format
        Does not perform checks so be carefull when used.
        
        Inputs:
        - fasta, full path to fasta file
        - outFile, full path to phylip output file
        
        Outputs:
        - outFile, fasta converted to phylip format
    '''
    
    # read fasta file
    names = []
    seq = []
    
    with open(outFile, 'w') as output:
        with open(fasta, 'r') as fasta_input:
            for line in fasta_input:
                if line.startswith('>'):
                    names.append(line.strip().replace('>',''))
                else:
                    seq.append(line.strip())
    
        # convert to phylip format
        # save to phylip file
        output.write('{} {}\n'.format(len(names), len(seq[0])))
        for i, d in enumerate(names):
            print d
            output.write('{}\n'.format(d))
            output.write('     {}\n\n'.format(seq[i]))
    
    
    
