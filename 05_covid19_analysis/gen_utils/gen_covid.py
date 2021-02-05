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
        ## B117 variants
    b117_variants = {
    'C3267T':'ORF1ab-T1001I',
    'C5388A':'ORF1ab-A1708D',
    'T6954C': 'ORF1ab-I2230T',
    'A23063T':'spike-N501Y',
    'C23271A':'spike-A570D',
    'C23604A':'spike-P681H',
    'C23709T':'spike-T716I',
    'T24506G':'spike-S982A',
    'G24914C':'spike-D1118H',
    'C27972T':'Orf8-Q27stop',
    'G28048T':'Orf8-R52I',
    'A28111G':'Orf8-Y73C',
    'G28280C':'N-D3L-80',
    'A28281T':'N-D3L-81',
    'T28282A':'N-D3L-82',
    'C28977T':'N-S235F'
    }
    
    ### b117 deletions
    b117_deletions ={
    "b117_del_11288": [11288,11289,11290,11291,11292,11293,11294,11295,11296],
    "b117_del_21765": [21765,21766,21767,21768,21769,21770],
    "b117_del_21991": [21991,21992,21993]
    }

    cal20_variants ={
    'A12878G':'ORF1ab-I4205V',
    'G17014T':'ORF1ab-D5584Y',
    'G21600T':'Spike-S13I',
    'G22018T':'Spike-W152C',
    'T22917G':'Spike-L452R',
    'A23403G':'Spike-D614G'
    }

    # south africa
    b1351_variants = {
    'C26456T':'E-P71L',
    'C28887T':'N-T205I',
    'G5230T':'ORF1a-K1655N',
    'A21801C':'Spike-D80A', 
    'C21614T':'Spike-L18F',
    'A22206G':'Spike-D215G',
    'G22813T':'Spike-K417N', 
    'G22299T':'Spike-R246I',
    'G23012A':'Spike-E484K',
    'C23664T':'Spike-A701V',
    'A23063T':'Spike-N501Y',
    'A23403G':'Spike-D614G'
    }

   # brazil p1
    b1128_variants = {
    'C3828T':'ORF1ab-S1188L',
    'A5648C':'ORF1ab-K1795Q',
    'G25912T':'ORF3a-G174C',
    'G28167A':'ORF8-E92K',
    'C28512G':'N-P80R',
    'C21614T':'Spike-L18F',
    'C21621A':'Spike-T20N',
    'C21638T':'Spike-P26S',
    'G21974T':'Spike-D138Y',
    'G22132T':'Spike-R190S',
    'G22813T':'Spike-K417N',
    'G23012A':'Spike-E484K',
    'A23063T':'Spike-N501Y',
    'C23525T':'Spike-H655Y',
    'C24642T':'Spike-T1027I'
    }	

    # brazil p2
    p2_variants = {
    'C3828T':'ORF1ab-S1188L',
    'A5648C':'ORF1ab-K1795Q',
    'G25912T':'ORF3a-G174C',
    'G28167A':'ORF8-E92K',
    'C28512G':'N-P80R',
    'C21614T':'Spike-L18F',
    'C21621A':'Spike-T20N',
    'C21638T':'Spike-P26S',
    'G21974T':'Spike-D138Y',
    'G22132T':'Spike-R190S',
    'G22813T':'Spike-K417N',
    'G23012A':'Spike-E484K',
    'A23063T':'Spike-N501Y',
    'C23525T':'Spike-H655Y',
    'C24642T':'Spike-T1027I',
    'C100T':'5UTR',
    'C28253T':'ORF8-F120F',
    'G28628T':'N-A119S', 
    'G28975T':'N-M234I',
    'C29754T':'3UTR'
    }
