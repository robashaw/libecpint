from lxml import etree

angular_shells = ["s", "p", "d", "f", "g", "h", "i"]

class Shell:
    def __init__(self, lval = 0, nexp = 0, nbfs = 0):
        self.lval = lval
        self.powers = []
        self.nexp = nexp
        self.exps = []
        self.contr = []

class Atom:
    def __init__(self, name="X", nshells = 0):
        self.name = name
        self.ncore = 0
        self.maxl = 0
        self.nshells = nshells
        self.shells = []
        
def parse_ecp(file):
    atoms = []
    atomnames = {}
    lines = file.readlines()
    nlines = len(lines)
    linenumber = 0
    atomnumber = 0
    
    while linenumber < nlines:
        line = lines[linenumber].strip()
        tokens = line.split(',')
        for token in tokens: 
            token = token.replace(' ', '')
            
        if (len(tokens) > 2):
            if (tokens[0].lower() == "ecp"):
                atom_name = tokens[1].lower()
                new_atom = Atom(name=atom_name)

                new_atom.ncore = int(tokens[2])
                new_atom.maxl = int(tokens[3])
                
                for i in range(new_atom.maxl+1):
                    linenumber += 1
                    line = lines[linenumber].strip()
                    
                    tokens = line.split(';')
                    for token in tokens:
                        token = token.replace(' ', '')
                        
                    if(len(tokens) > 1):
                        nprims = int(tokens[0])
                        l = i-1
                        if (i==0):
                            l = new_atom.maxl
                        new_shell = Shell(lval=l, nexp=nprims, nbfs=1)
                        
                        for token in tokens[1:]:
                            subtokens = token.split(',')
                            if (len(subtokens) == 3):
                                new_shell.powers.append(subtokens[0])
                                new_shell.exps.append(subtokens[1])
                                new_shell.contr.append(subtokens[2])
                        
                        new_atom.shells.append(new_shell)
                
                atoms.append(new_atom)
                
        linenumber += 1
    
    return atoms

def write_ecp_basis(atoms, name):
    filename = "xml/" + name + ".xml"
    with open(filename, 'wb') as new_file:
        root = etree.Element("root", name=name)
        tree = etree.ElementTree(root)
        
        for atom in atoms:
            child = etree.SubElement(root, atom.name, ncore = str(atom.ncore), maxl=str(atom.maxl))
            
            for shell in atom.shells:
                schild = etree.SubElement(child, "Shell", lval=str(shell.lval), nexp=str(shell.nexp))
                
                for i in range(shell.nexp):
                    try:
                        xchild = etree.SubElement(schild, "nxc", n=shell.powers[i], x=shell.exps[i], c=shell.contr[i])
                    except:
                        print("ERROR in " + filename + " (" + name + ") : atom " + atom.name + ", shell type " + shell.lval)
                        print("Expected no. of powers/exps/coeffs:", shell.nexp)
                        print("Actual no. of powers: ", len(shell.powers))
                        print("Actual no. of exps: ", len(shell.exps))
                        print("Actual no. of coeffs: ",  len(shell.contr))
        tree.write(new_file, pretty_print = True)

if __name__ == "__main__":
    import sys
    name = sys.argv[1]
    input_file = open('raw/' + name + '.ecp', 'r')
    atoms = parse_ecp(input_file)
    write_ecp_basis(atoms, name)
    
