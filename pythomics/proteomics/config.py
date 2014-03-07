#We use X!Tandem notation for enzymatic digestion:
#http://thegpm.org/TANDEM/api/pcs.html
ENZYMES = {'trypsin': '[KR]|{P}', 'lysc': '[K]|[X]', 'none': '|',
           'argc': '[R]|[X]', 'gluc': '[E]|[X]', 'aspn': '[X]|[D]',
           'V8': '[ED]|[X]'}

#masses from: http://www.weddslist.com/ms/tables.html#tm4
#stored as mass, charge
RESIDUE_MASSES =  {'G': (57.021464, 0),
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