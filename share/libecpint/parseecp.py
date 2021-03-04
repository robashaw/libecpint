# Script to parse ECP files into XML format for use by libecpint
# Usage: python parseecp.py [NAME]
# will convert raw/NAME.ecp into xml/NAME.xml
# which can then be found by libecpint using reference "NAME"
from lxml import etree

# letter symbols for angular momentum shells
# can be added to if needed
angular_shells = ["s", "p", "d", "f", "g", "h", "i"]


class Shell:
    """ Container for a shell of an ECP
        i.e. the Gaussian expansion of a fixed angular momentum
        
        Members:
        lval - angular momentum quantumm number of shell
        powers - the power of r for each Gaussian in the shell
        nexp - the number of primitive Gaussians (exponents) in the shell
        exps - the exponents of the Gaussians
        contr - the contraction coefficients of the Gaussians
    """
    def __init__(self, lval = 0, nexp = 0):
        """Creates a  Shell with angular momentum lval and nexp exponents"""
        self.lval = lval
        self.powers = []
        self.nexp = nexp
        self.exps = []
        self.contr = []

class Atom:
    """ Container for an ECP for a specific atom type
        
        Members:
        name - name of the atom, e.g. 'O' for oxygen
        ncore - the number of core electrons in the ECP
        maxl - the maximum angular momentum shell in the ECP
        nshells - the number of shells in the ECP
        shells - an array of Shell objects describing the shells in the ECP
    """
    def __init__(self, name="X", nshells = 0):
        """Creates a blank Atom with the given name"""
        self.name = name
        self.ncore = 0
        self.maxl = 0
        self.nshells = nshells
        self.shells = []
        
def tokenize(line, sep=','):
    """Given a line of input, cleans up and returns tokens,
       split by the separator (sep)
    """
    # strip out whitespace and split by separator
    line = line.strip()
    tokens = line.split(sep)
    # get rid of additional whitespace
    for token in tokens: 
        token = token.replace(' ', '')
    return tokens
        
def parse_ecp(file):
    """Given a MOLPRO-format ECP file, returns Atom objects (ECPs) for every
       atom type defined in that file.  
    """
    atoms = []
    atomnames = {}
    
    # read in the file
    lines = file.readlines()
    nlines = len(lines)
    
    # loop to the end of the file
    linenumber = 0
    atomnumber = 0
    while linenumber < nlines:
        tokens = tokenize(lines[linenumber], sep=',')
        
        # all lines where ECP definitions start will have 
        # at least two bits "ecp,AtomName"
        if (len(tokens) > 2):
            if (tokens[0].lower() == "ecp"):
                # Found an ECP definition
                atom_name = tokens[1].lower()
                # Create a container for the ECP
                new_atom = Atom(name=atom_name)

                # the next two tokens should be the no. of core electrons
                # and the maximum angular momentum shell of the ECP
                new_atom.ncore = int(tokens[2])
                new_atom.maxl = int(tokens[3])
                # there should then be maxl+1 lines definining the shells
                # in the order [maxL, 0, 1, ..., maxL-1]                
                for i in range(new_atom.maxl+1):
                    linenumber += 1
                    tokens = tokenize(lines[linenumber], sep=';')
                    
                    # could be blank lines or comments, so check first
                    if(len(tokens) > 1):
                        # shell definition has form "nx; n,x,c; n,x,c; ..."
                        # defining nx Gaussians of the form c * r^n * exp(-x*r^2)
                        nprims = int(tokens[0])
                        l = i-1
                        if (i==0):
                            l = new_atom.maxl
                        
                        # create a container for the shell
                        new_shell = Shell(lval=l, nexp=nprims)
                        # fill in the details of the shell as described above
                        for token in tokens[1:]:
                            subtokens = token.split(',')
                            if (len(subtokens) == 3):
                                new_shell.powers.append(subtokens[0])
                                new_shell.exps.append(subtokens[1])
                                new_shell.contr.append(subtokens[2])
                        
                        # append the new shell to the Atom
                        new_atom.shells.append(new_shell)
                # end of Atom definition, append
                atoms.append(new_atom)
                
        linenumber += 1
    
    # Return all the atoms found
    return atoms

def write_ecp_basis(atoms, name):
    """Given a list of Atom objects defining ECPs, and a name for the ECP basis,
       this writes the basis to XML file. 
    """
    filename = "xml/" + name + ".xml"
    
    # write into a binary xml file using LXML package
    with open(filename, 'wb') as new_file:
        # format is Root -> Atom1 --> Shell1
        #                            Shell2 ... etc
        #                   Atom2 --> etc.
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
                        # something is wrong with the definition of this shell
                        print("ERROR in " + filename + " (" + name + ") : atom " + atom.name + ", shell type " + shell.lval)
                        print("Expected no. of powers/exps/coeffs:", shell.nexp)
                        print("Actual no. of powers: ", len(shell.powers))
                        print("Actual no. of exps: ", len(shell.exps))
                        print("Actual no. of coeffs: ",  len(shell.contr))
        tree.write(new_file, pretty_print = True)

if __name__ == "__main__":
    # Given a raw MOLPRO ECP file in raw/name.ecp, writes xml/name.xml
    import sys
    name = sys.argv[1]
    input_file = open('raw/' + name + '.ecp', 'r')
    atoms = parse_ecp(input_file)
    write_ecp_basis(atoms, name)
    
