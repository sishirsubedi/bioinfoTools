import os
import sys

class covid_basics:

    reference = "MN908947"

    AA_TriToS={
        "Ala": "A",
        "Asx": "B",
        "Cys": "C",
        "Asp": "D",
        "Glu": "E",
        "Phe": "F",
        "Gly": "G",
        "His": "H",
        "Ile": "I",
        "Xle": "J",
        "Lys": "K",
        "Leu": "L",
        "Met": "M",
        "Asn": "N",
        "Hyp": "O",
        "Pro": "P",
        "Gln": "Q",
        "Arg": "R",
        "Ser": "S",
        "Thr": "T",
        "Glp": "U",
        "Val": "V",
        "Trp": "W",
        "Ter": "X",
        "Tyr": "Y",
        "Glx": "Z"}

    nsp12_region_start = 13442
    nsp12_region_stop = 16237
    s_region_start = 21563
    s_region_stop = 25385

class covid_variants:
    
    snvs_threshold ={
    "B_1_1_7" : 10,
    "B_1_429" : 10,
    "B_1_427" : 9,
    "B_1_351" : 8,
    "P_1" : 15,
    "P_2" : 10,
    "P_3" : 13,
    "B_1_525" : 10,
    "B_1_526" : 10,
    "B_1_2" : 9
    }
    
def get_covid_variants():
    
    var_dict =  {

    'B_1_1_7' : {
    'C3267T':'ORF1ab-T1001I',
    'C5388A':'ORF1ab-A1708D',
    'T6954C': 'ORF1ab-I2230T',
    'A23063T':'S-N501Y',
    'C23271A':'S-A570D',
    'A23403G':'S-D614G',
    'C23604A':'S-P681H',
    'C23709T':'S-T716I',
    'T24506G':'S-S982A',
    'G24914C':'S-D1118H',
    'C27972T':'ORF8-Q27*',
    'G28048T':'ORF8-R52I',
    'A28111G':'ORF8-Y73C',
    'G28280C':'N-D3H',
    'A28281T':'N-D3V',
    'T28282A':'N-D3E',
    'C28977T':'N-S235F'
    },

    'B_1_429' : {
    'C1059T':'ORF1ab-T265I', 
    'A12878G':'ORF1ab-I4205V', 
    'C14408T':'ORF1ab-P4715L', 
    'G17014T':'ORF1ab-D5584Y',  
    'G21600T':'S-S13I', 
    'G22018T':'S-W152C', 
    'T22917G':'S-L452R', 
    'A23403G':'S-D614G',  
    'G25563T':'ORF3a-Q57H',  
    'C28887T':'N-T205I'
    },

    'B_1_427' : {
    'C1059T':'ORF1ab-T265I', 
    'G9738C':'ORF1ab-S3158T',  
    'C16394T':'ORF1ab-P5377L', 
    'G17014T':'ORF1ab-D5584Y',
    'G21600T':'S-S13I',  
    'T22917G':'S-L452R', 
    'A23403G':'S-D614G', 
    'G25563T':'ORF3a-Q57H',  
    'C28887T':'N-T205I'
    },

    'B_1_351'  : {
    'C26456T':'E-P71L',
    'C28887T':'N-T205I',
    'G5230T':'ORF1ab-K1655N',
    'G22813T':'S-K417N', 
    'G23012A':'S-E484K',
    'A23063T':'S-N501Y',
    'A23403G':'S-D614G',
    'C23664T':'S-A701V'
    },

    'P_1' : {
    'C3828T':'ORF1ab-S1188L', 
    'A5648C':'ORF1ab-K1795Q', 
     'C21614T':'S-L18F', 
     'C21621A':'S-T20N', 
     'C21638T':'S-P26S', 
     'G21974T':'S-D138Y', 
     'G22132T':'S-R190S', 
     'A22812C':'S-K417T', 
     'G23012A':'S-E484K', 
     'A23063T':'S-N501Y', 
     'A23403G':'S-D614G', 
     'C23525T':'S-H655Y', 
     'C24642T':'S-T1027I', 
     'G28167A':'ORF8-E92K', 
     'C28512G':'N-P80R', 
    },
    # MISSING--
    # ORF1ab	I760T
    # ORF9	Q77E
    # ORF3a	C174G
    # ORF1ab	F681L
    # ORF1ab	E5662D

    'P_2' : {
    'T10667G':'ORF1ab-L3468V', 
    'C12053T':'ORF1ab-L3930F',  
    'C14408T':'ORF1ab-P4715L', 
    'G23012A':'S-E484K',  
    'A23403G':'S-D614G', 
    'G25088T':'S-V1176F',  
    'G28628T':'N-A119S', 
    'G28881A':'N-R203K', 
    'G28883C':'N-G204R', 
    'G28975T':'N-M234I'
    },

    'P_3' : {
    'A4962G':'ORF1ab-D1554G', 
    'T9867C':'ORF1ab-L3201P',  
    'C11308A':'ORF1ab-D3681E', 
    'C12053T':'ORF1ab-L3930F', 
    'C17339T':'ORF1ab-A5692V', 
    'G23012A':'S-E484K',  
    'A23063T':'S-N501Y', 
    'A23403G':'S-D614G',  
    'C23604A':'S-P681H',  
    'G24836A':'S-E1092K',  
    'C24863T':'S-H1101Y',  
    'G25088T':'S-V1176F',
    'A27897C':'ORF8-K2Q'
    },

    'B_1_525' : {
    'C6285T':'ORF1ab-T2007I',  
    'C14408T':'ORF1ab-P4715L', 
    'C21762T':'S-A67V', 
    'G23012A':'S-E484K', 
    'A23403G':'S-D614G', 
    'G23593C':'S-Q677H', 
    'T24224C':'S-F888L',
    'T26767C':'M-I82T',
    'C28308G':'N-A12G', 
    'C28887T':'N-T205I'
    },

    #https://www.medrxiv.org/content/10.1101/2020.11.30.20241265v2.full-text
    # B_1_177_variants = { 
    #     'C14408T':'ORF1ab-P4715L', 
    #     'C22227T':'S-A222V', 
    #     'A23403G':'S-D614G' ,
    #     'C28932T':'N-A220V',
    #     'G29645T': 'ORF10-V30L'
    # }

    # B_1_1_64_variants = {
    #     'C14408T':'ORF1ab-P4715L', 
    #     'G21724T':'S-L54F' ,
    #     'A23403G':'S-D614G' ,
    #     'G28881A':'N-R203K',
    #     'G28883A' :'N-G204R'
    # }

    'B_1_526' : {
    'C1059T':'ORF1ab-T265I',
    'T9867C':'ORF1ab-L3201P',
    'C14408T': 'ORF1ab-P4715L',
    'A16500C':'ORF1ab-Q5412H',
    'C21846T':'S-T95I',
    'A22320G':'S-D253G',
    'A23403G':'S-D614G',
    'G25563T':'ORF3a-Q57H',
    'C25517T':'ORF3a-P42L',
    'C27925T':'ORF8-T11I'
    },
    
    'B_1_2' : {
    'C1059T':'ORF1ab-T265I',
    'C10319T':'ORF1ab-L3352F',
    'C14408T': 'ORF1ab-P4715L',
    'A23403G':'S-D614G',
    'A18424G':'ORF1ab-N6054D',
    'C21304T':'ORF1ab-R7014C',
    'G25907T':'ORF3a-G172V',
    'C28472T':'N-P67S',
    'C28869T':'N-P199L'
    }
    }

    return var_dict

def get_s_domain(loc):

    domain = ""

    protein_domain={
        'S1-NTD':'16-305',
        'S1-RBD':'330-521',
        'Furin Cleavage':'682-685',
        'S2-FP':'816-833',
        'S2-HR1':'908-985',
        'S2-CH':'986-1035',
        'S2-CD':'1076-1141'}

    if loc<=681:
        domain = 'S1'
    elif loc>=686:
        domain = 'S2'

    for d in protein_domain:
        start = int(protein_domain[d].split('-')[0])
        end = int(protein_domain[d].split('-')[1])

        if loc>=start and loc<=end:
            domain = d
    
    return domain