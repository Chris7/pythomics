HYDROGEN = 1.007825
PROTON = 1.00727
#NEUTRON = 1.008701
# we use the mass of a neutron of c13
NEUTRON = 1.003355
# Neutrons with mass defects for a given element
CARBON_NEUTRON = 1.003355
HYDROGEN_NEUTRON = 1.006277
NITROGEN_NEUTRON = 0.997035

#We use X!Tandem notation for enzymatic digestion:
#http://thegpm.org/TANDEM/api/pcs.html
ENZYMES = {'trypsin': '[KR]|{P}', 'lysc': '[K]|[X]', 'none': '|',
           'argc': '[R]|[X]', 'gluc': '[E]|[X]', 'aspn': '[X]|[D]',
           'V8': '[ED]|[X]'}

#masses from: http://www.weddslist.com/ms/tables.html#tm4
#stored as mass, charge. Note these are the masses after a water loss
RESIDUE_MASSES =  {
  'G': (57.021464, 0),
  'A': (71.037114, 0),
  'S': (87.032028, 0),
  'P': (97.052764, 0),
  'V': (99.068414, 0),
  'T': (101.047678, 0),
  'C': (103.009184, 0),
  'I': (113.084064, 0),
  'L': (113.084064, 0),
  'N': (114.042927, 0),
  'D': (115.026943, 0),
  'Q': (128.058578, 0),
  'K': (128.094963, 1),
  'E': (129.042593, 0),
  'M': (131.040485, 0),
  'H': (137.058912, 1),
  'F': (147.068414, 0),
  'R': (156.101111, 1),
  'Y': (163.063329, 0),
  'W': (186.079313, 0)
}

RESIDUE_COMPOSITION = {
    'G': {'C': 2, 'H': 5, 'O': 2, 'N': 1},
    'A': {'C': 3, 'H': 7, 'O': 2, 'N': 1},
    'S': {'C': 3, 'H': 7, 'O': 3, 'N': 1},
    'P': {'C': 5, 'H': 9, 'O': 2, 'N': 1},
    'V': {'C': 5, 'H': 11, 'O': 2, 'N': 1},
    'T': {'C': 4, 'H': 9, 'O': 3, 'N': 1},
    'C': {'C': 3, 'H': 7, 'O': 2, 'N': 1, 'S': 1},
    'I': {'C': 6, 'H': 13, 'O': 2, 'N': 1},
    'L': {'C': 6, 'H': 13, 'O': 2, 'N': 1},
    'N': {'C': 4, 'H': 8, 'O': 3, 'N': 2},
    'D': {'C': 4, 'H': 7, 'O': 4, 'N': 1},
    'Q': {'C': 5, 'H': 10, 'O': 3, 'N': 2},
    'K': {'C': 6, 'H': 14, 'O': 2, 'N': 2},
    'E': {'C': 5, 'H': 9, 'O': 4, 'N': 1},
    'M': {'C': 5, 'H': 11, 'O': 2, 'N': 1, 'S': 1},
    'H': {'C': 6, 'H': 9, 'O': 2, 'N': 3},
    'F': {'C': 9, 'H': 11, 'O': 2, 'N': 1},
    'R': {'C': 6, 'H': 14, 'O': 2, 'N': 4},
    'Y': {'C': 9, 'H': 11, 'O': 3, 'N': 1},
    'W': {'C': 11, 'H': 12, 'O': 2, 'N': 2},
}

ELEMENTS = {
    'C': 12.0,
    'H': HYDROGEN,
    'N': 14.003074,
    'O': 15.994915,
    'S': 31.9720707,
}

PI_RESIDUE_CHARGE = {'D': (3.65, -1.0),
                     'E': (4.25, -1.0),
                     'C': (8.18, -1.0),
                     'Y': (10.07, -1.0),
                     'H': (6.00, 1.0),
                     'K': (10.53, 1.0),
                     'R': (12.48, 1.0)}

#stored as mass, charge
MODIFICATION_MASSES = {'h': (1.007825, 0),
               'h2o': (18.010565, 0),
               'h2o': (18.0106, 0),
               'nh3': (17.026549, 0),
               'ch2': (14.015650, 0),
               'methylation': (14.015650, 0),
               'o': (15.994915, 0),
               'oxidation': (15.994915, 0),
               'acetylation': (42.010565, 0),
               'acetylation': (42.0106, 0),
               'carbamidation': (57.021464, 0),
               'carboxylation': (58.005479, 0),
               'phosphorylation': (79.966330, 0),
               'amidation': (0.984016, 0),
               'formylation': (27.994915, 0),
               'cho': (29.002739665, 0),
               'nh2': (16.01872407, 0),
               'co': (27.99491463, 0),
               'oh': (17.00274, 0)
               }

#neutral losses, shown as residue, mass loss, fragment type affected, and display name
LOSS_MASSES = {'K': ((-1*MODIFICATION_MASSES['nh3'][0],('a','b','y'),'-nh3'),),
              'R': ((-1*MODIFICATION_MASSES['nh3'][0],('a','b','y'),'-nh3'),),
              'Q': ((-1*MODIFICATION_MASSES['nh3'][0],('a','b','y'),'-nh3'),),
              'N': ((-1*MODIFICATION_MASSES['nh3'][0],('a','b','y'),'-nh3'),),
              'S': ((-1*MODIFICATION_MASSES['h2o'][0],('b'),'-h2o'),),
              'T': ((-1*MODIFICATION_MASSES['h2o'][0],('b'),'-h2o'),),
              'D': ((-1*MODIFICATION_MASSES['h2o'][0],('b','y'),'-h2o'),),
              'E': ((-1*MODIFICATION_MASSES['h2o'][0],('b','y'),'-h2o'),)}

LABEL_SCHEMES = {
    'O18': {
        'O18': {4.008492: set([']'])},
        'O16': {0: set([])},
    },
    'N15': {
        'N15': {0.997035: set(['X']), 1.99407: set(['K', 'Q', 'N']), 2.991105: set(['H']), 3.98814: set(['R'])},
        'N14': {0: set([])},
        },
    'mTRAQ4':{
        'Medium': {4.0071: set([']', 'K']), 4.0071*2: set(['R'])},
        'Light': {0: set([])}
    },
    'mTRAQ48':{
        'Heavy': {8.014199: set([']', 'K']), 8.014199*2: set(['R'])},
        'Medium': {4.0071: set([']', 'K']), 4.0071*2: set(['R'])},
        'Light': {0: set([])}
    },
    'mTRAQ8':{
        'Heavy': {8.014199: set([']', 'K']), 8.014199*2: set(['R'])},
        'Light': {0: set([])}
    },
    'K6R6': {
              'Light': {0: set([])},
              'Heavy': {6.02013: set(['R']), 4.02511: set(['K'])}
             },
    'K8R10': {'Heavy': {8.0142: set(['K']), 10.00827: set(['R'])},
              'Light': {0: set([])}
              },
    'K4K8R6R10': {'Heavy': {8.0142: set(['K']), 10.00827: set(['R'])},
                  'Light': {0: set([])},
                  'Medium': {6.02013: set(['R']), 4.02511: set(['K'])}
                  },
    'TMT6': {
        '126': {126.127725: None},
        '127': {127.124760: None},
        '128': {128.134433: None},
        '129': {129.131468: None},
        '130': {130.141141: None},
        '131': {131.138176: None},
    },
    'TMT10': {
        '126': {126.127725: None},
        '127N': {127.124760: None},
        '127C': {127.131079: None},
        '128N': {128.128114: None},
        '128C': {128.134433: None},
        '129N': {129.131468: None},
        '129C': {129.137787: None},
        '130N': {130.134822: None},
        '130C': {130.141141: None},
        '131': {131.138176: None},
    },
    'iTRAQ8': {
        '117': {117.1: None},
        '116': {116.1: None},
        '115': {115.1: None},
        '114': {114.1: None},
        '113': {113.1: None},
        '118': {118.1: None},
        '119': {119.1: None},
        '121': {121.1: None},
    },
    'iTRAQ4': {
        '117': {117.1: None},
        '116': {116.1: None},
        '115': {115.1: None},
        '114': {114.1: None},
    }
}