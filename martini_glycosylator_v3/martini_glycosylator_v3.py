import sys
import pandas as pd
import glob, os
import csv
import argparse

itp_files = glob.glob('gly_atoms*.itp') + glob.glob('gly_bonds*.itp') + glob.glob('gly_angles*.itp') + glob.glob('gly_dihedrals*.itp') + glob.glob('gly_constaints*.itp')

for f in glob.glob("*_gly.pdb"):
    os.remove(f)
for f in glob.glob("*chain*pdb"):
    os.remove(f)
for f in itp_files:
    os.remove(f)

parser = argparse.ArgumentParser(description='***********Martini glycosylator*************')

parser.add_argument("-d", default=1, type=float, dest='dist', help="The glycan is placed  d*3.5 angstrom away from asparagine. A value between 1 to 2 works for most cases")
parser.add_argument("-c", required=True, type=str, dest='pdb', help="CG Protein PDB file with chains specified is required")
parser.add_argument("-f", required=True, type=str, dest='dat', help="Input file with ASN_site, chain and glycan attached is required")

args = parser.parse_args()
dist = args.dist
print('d = {}\n'.format(dist))
protein = args.pdb
print('protein = {}\n'.format(protein))
input_file = args.dat
print('input data: \n')


fileDir = os.path.dirname(os.path.realpath('__file__'))
input_dat = pd.read_csv(input_file, delim_whitespace=True)
print(input_dat.to_string(index=False)+'\n')
input_dat_temp = pd.read_csv(input_file, delim_whitespace=True)
input_pdb = fileDir+ '/all_gly_pdbs/'+  input_dat['glycan'] + '.pdb'
input_dat['Prot_PDB'] = ''
input_dat['glycan_itp'] = fileDir+ '/all_gly_itps/'+  input_dat['glycan'] + '.itp'
input_dat['glycan_pdb'] = fileDir+ '/all_gly_pdbs/'+  input_dat['glycan'] + '.pdb'
input_dat['BB'] = ''
input_dat['last_index'] = ''
input_dat['last_resid'] = ''
protein_path = fileDir+'/'+protein

################# for ITP side of things ####################
############################################################# 
#count number of lines in a file
def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

## splits the protein PDB into chains ###
def chain_splitter(chains, chain_info, pdb, protein, input_dat):
    for i in (unique_chains):
        count  = int(0)
        with open(protein_path, 'r') as pdb, open(protein[:-4] + '_chain_'+ i +'.pdb', 'a') as outfile, open(protein[:-4] + '_chain_'+ i +'_gly.pdb', 'a') as glyoutfile:
            clean_words = ['MODEL', 'REMARK', 'TITLE', 'CRYST1', 'TER', 'ENDMDL', 'END']
            for line in pdb:
                if not any(clean_word in line for clean_word in clean_words):
                    pdb_line = line.split()
                    if pdb_line[4][:1] == i:
                        outfile.write(line)
                        glyoutfile.write(line)
                        count = count + 1
                        last_resid = pdb_line[5]
			if len(pdb_line[4]) > 1: #if resid goes beyond 999 then A999, resid is column - first character
		           #print(pdb_line) 
			   last_resid = pdb_line[4][1:]
        print('Number of atoms in chain {} = {}'.format(i, count))
        chain_info = chain_info.append({'chain': i , 'atoms': count},  ignore_index=True)
 
        for j in range(0, input_dat.shape[0]):
            if input_dat.loc[j,'chain'] == i:
                
                input_dat.loc[j,'last_index']= int(count)
                input_dat.loc[j, 'last_resid'] = int(last_resid)
                
    #update input_dat file
    for i in range(0, input_dat.shape[0]):
        input_dat.loc[i,'Prot_PDB'] = protein[:-4]+ '_chain_'+ input_dat.loc[i,'chain']+'.pdb'
#         input_dat.loc[i, 'last_resid'] =
    return input_dat, chain_info
#       print('wrote {}'.format(glyoutfile.name))

#reads protein pdb file and returns BB atom indexes and last resid and index for given glycosylation sites## 
def protpdb_reader(glycan_data):
    print('\nFinding backbone indexes for asparagine sites')
    clean_words = ['MODEL', 'REMARK', 'TITLE', 'CRYST1', 'TER', 'ENDMDL','END']
    last_resid = 0
    last_index = 0
#     glycan_data['BB'] = ''
    for i in range (0, glycan_data.shape[0]):
        with open (glycan_data.loc[i,'Prot_PDB']) as pdb_data:
            for line in pdb_data:
                pdb_line = line.split()
                if not any(clean_word in line for clean_word in clean_words):
                    #if resid is more than 999
                    #print(pdb_line[4])
                    if (len(pdb_line[4]) > 1 and int(pdb_line[4][1:])== glycan_data.loc[i,'ASN_site'] and pdb_line[4][:1]== glycan_data.loc[i,'chain'] and pdb_line[2] == 'BB'):
                        print("Backbone index for chain {} ASN {} is {}".format(glycan_data.loc[i,'chain'],glycan_data.iloc[i,0], pdb_line[1]))
			glycan_data.loc[i, 'BB'] = int(pdb_line[1])
			chain = glycan_data.loc[i,'chain']
                       	break

                    elif (len(pdb_line[4]) == 1 and int(pdb_line[5])== glycan_data.loc[i,'ASN_site'] and pdb_line[4]== glycan_data.loc[i,'chain'] and pdb_line[2] == 'BB'):
                        print("Backbone index for chain {} ASN {} is {}".format(glycan_data.loc[i,'chain'],glycan_data.iloc[i,0], pdb_line[1]))
			glycan_data.loc[i, 'BB'] = int(pdb_line[1])
			chain = glycan_data.loc[i,'chain']
                        break
    return(glycan_data) 
######## renumbers the BB atoms
def set_bb(input_dat, chain_info):
	for i in range(0,chain_info.shape[0]):
    		correction = 0
    		if i == 0:
        		correction = 0
			chain_info.loc[i, 'correction'] = int(correction)
    		else:
        		for j in range(0, i):
            			correction = correction + int(chain_info.loc[j,'atoms'])
        		chain_info.loc[i, 'correction'] = int(correction)

	for i in range(0,input_dat.shape[0]):
		for j in range(0,chain_info.shape[0]):
			if input_dat.loc[i,'chain']== chain_info.loc[j,'chain']:
				input_dat.loc[i,'BB'] = input_dat.loc[i,'BB'] - chain_info.loc[j,'correction']
	#print(input_dat)
	return input_dat

#reads glycan itp files
def itp_reader(filename):
    data_file = filename
    data_file_delimiter = ' '
    column_names = [i for i in range(0, 8)]
    df = pd.read_csv(data_file, header=None, delim_whitespace=True, names=column_names)
    df = df.fillna('')
    return df

#Reads ITP files, corrects indexes and adds Asparagine connections####################################
def itp_indexer(df, last_pdb_index, last_resid, asn_par, chain):
    add_index = 0
    for i in range(0, df.shape[0]):
        if df.iloc[i,0] == '[atoms]':
            j = i+1
            while df.iloc[j,0] != '[bonds]':
#                 print(df.iloc[j,0])
                df.iloc[j,0] = int(df.iloc[j,0]) + last_pdb_index
                df.iloc[j,2] = int(df.iloc[j,2]) + last_resid 
                new_last_index = df.iloc[j,0]
                new_last_resid = df.iloc[j,2]
                df.loc[[j]].to_csv('gly_atoms_'+ chain + '.itp', sep='\t', index = False, header = False, mode='a')
                j = j + 1
                add_index = add_index + 1

        if df.iloc[i,0] == '[bonds]':
            j = i+1
            pd.DataFrame([asn_par['bond1']]).to_csv('gly_bonds_'+ chain + '.itp', sep='\t', index = False, header = False, mode='a')
            while df.iloc[j,0] not in ['[constraints]','[angles]']:
                df.iloc[j,0] = int(df.iloc[j,0]) + last_pdb_index
                df.iloc[j,1] = int(df.iloc[j,1]) + last_pdb_index
                df.iloc[j,2] = int(df.iloc[j,2]) 
                df.loc[[j]].to_csv('gly_bonds_'+ chain + '.itp', sep='\t', index = False, header = False, mode='a')
                j = j + 1 

        if df.iloc[i,0] == '[constraints]':
            j = i+1
            while df.iloc[j,0] != '[angles]':
                df.iloc[j,0] = int(df.iloc[j,0]) + last_pdb_index
                df.iloc[j,1] = int(df.iloc[j,1]) + last_pdb_index
		df.iloc[j,2] = int(df.iloc[j,2])
                df.loc[[j]].to_csv('gly_constraints_'+ chain + '.itp', sep='\t', index = False, header = False, mode='a')
                j = j + 1 

        if df.iloc[i,0] == '[angles]':
            pd.DataFrame([asn_par['angle1']]).to_csv('gly_angles_'+ chain + '.itp', sep='\t', index = False, header = False, mode='a')
            pd.DataFrame([asn_par['angle2']]).to_csv('gly_angles_'+ chain + '.itp', sep='\t', index = False, header = False, mode='a')
            pd.DataFrame([asn_par['angle3']]).to_csv('gly_angles_'+ chain + '.itp', sep='\t', index = False, header = False, mode='a')
            j = i+1
            while df.iloc[j,0] != '[dihedrals]':
                df.iloc[j,0] = int(df.iloc[j,0]) + last_pdb_index
                df.iloc[j,1] = int(df.iloc[j,1]) + last_pdb_index
                df.iloc[j,2] = int(df.iloc[j,2]) + last_pdb_index
                df.iloc[j,5] = int(df.iloc[j,5]) 
                df.loc[[j]].to_csv('gly_angles_'+ chain + '.itp', sep='\t', index = False, header = False, mode='a')
                j = j + 1 

        if df.iloc[i,0] == '[dihedrals]':
            pd.DataFrame([asn_par['dihedral1']]).to_csv('gly_dihedrals_'+ chain + '.itp', sep='\t', index = False, header = False, mode='a')
            pd.DataFrame([asn_par['dihedral2']]).to_csv('gly_dihedrals_'+ chain + '.itp', sep='\t', index = False, header = False, mode='a')
            j = i+1
            while j != df.shape[0]:
                df.iloc[j,0] = int(df.iloc[j,0]) + last_pdb_index
                df.iloc[j,1] = int(df.iloc[j,1]) + last_pdb_index
                df.iloc[j,2] = int(df.iloc[j,2]) + last_pdb_index
                df.iloc[j,3] = int(df.iloc[j,3]) + last_pdb_index
                df.iloc[j,7] = int(df.iloc[j,7])
                df.loc[[j]].to_csv('gly_dihedrals_'+ chain + '.itp', sep='\t', index = False, header = False, mode='a')
                j = j + 1          
    return(new_last_index, new_last_resid)   

################# for PDB side of things ####################
############################################################# 
#reads protein pdb file and returns BB atom indexes for given glycosylation sites## 
def glypdb_reader(pdb, last_pdb_index, last_resid):
    clean_words = ['MODEL', 'REMARK', 'TITLE', 'CRYST1', 'TER', 'ENDMDL', 'END']
    with open (pdb) as oldfile, open('glycans.pdb', 'a') as newfile:
        for line in oldfile:
            if not any(clean_word in line for clean_word in clean_words):
                pdb_line = line.split()
                pdb_line[1] = last_pdb_index + 1
                last_pdb_index = last_pdb_index + 1
                pdb_line[4] = int(pdb_line[4]) + last_resid
                newfile.write("{:>4}{:>7}{:>5}{:>4}{:>6}{:>12}{:>8}{:>8}{:>6}{:>6}\n".format(*pdb_line))

##reads protein_pdb files and returns coordinates of BB and side chains of all ASN sites in input.dat #########                               
def protpdb_coord_reader(pdb, glycan_data, gly_num):
    coord = pd.DataFrame(columns=['site','x1','y1','z1','x2','y2','z2'])
    clean_words = ['MODEL', 'REMARK', 'TITLE', 'CRYST1', 'TER', 'ENDMDL','END']
    with open (pdb) as pdb_data:
        for line in pdb_data:
            pdb_line = line.split() 
            if not any(clean_word in line for clean_word in clean_words):
#                 print(line)
#                 for i in range(0, glycan_data.shape[0]):
                
                if (len(pdb_line[4])==1 and int(pdb_line[5])== glycan_data.loc[gly_num,'ASN_site'] and pdb_line[4]== glycan_data.loc[gly_num,'chain'] and pdb_line[2] == 'BB'):
                    print('Found atom {}-N{}-chain {}'.format(pdb_line[2],pdb_line[5],pdb_line[4]))
                    bb_coord = (pdb_line[6],pdb_line[7],pdb_line[8])
                #if resid mores than 999
                elif (len(pdb_line[4])>1 and int(pdb_line[4][1:])== glycan_data.loc[gly_num,'ASN_site'] and pdb_line[4][:1]== glycan_data.loc[gly_num,'chain'] and pdb_line[2] == 'BB'):
                    print('Found atom {}-N{}-chain {}'.format(pdb_line[2],pdb_line[5],pdb_line[4]))
                    bb_coord = (pdb_line[6],pdb_line[7],pdb_line[8])

                        
                elif (len(pdb_line[4])==1 and int(pdb_line[5])== glycan_data.loc[gly_num,'ASN_site'] and pdb_line[4]== glycan_data.loc[gly_num,'chain'] and pdb_line[2] == 'SC1'):
                    print('Found atom {}-N{}-chain {}'.format(pdb_line[2],pdb_line[5],pdb_line[4]))                    
                    sc1_coord = (pdb_line[6],pdb_line[7],pdb_line[8])
                    coord = coord.append({'site':glycan_data.iloc[gly_num,0],'x1': bb_coord[0], 'y1': bb_coord[1], 'z1': bb_coord[2], 'x2': sc1_coord[0], 'y2': sc1_coord[1], 'z2': sc1_coord[2]}, ignore_index=True)
                #if resid more than 999    
                elif (len(pdb_line[4])>1 and int(pdb_line[4][1:])== glycan_data.loc[gly_num,'ASN_site'] and pdb_line[4][:1]== glycan_data.loc[gly_num,'chain'] and pdb_line[2] == 'SC1'):
                    print('Found atom {}-N{}-chain {}'.format(pdb_line[2],pdb_line[5],pdb_line[4]))
                    sc1_coord = (pdb_line[6],pdb_line[7],pdb_line[8])
                    coord = coord.append({'site':glycan_data.iloc[gly_num,0],'x1': bb_coord[0], 'y1': bb_coord[1], 'z1': bb_coord[2], 'x2': sc1_coord[0], 'y2': sc1_coord[1], 'z2': sc1_coord[2]}, ignore_index=True)


    return(coord) 

##### finds ASN coordinates and moves glycan in the direction point BB>SC
def gly_to_prot_attch(input_dat, dist, input_dat_temp):
    for i in range (0, input_dat.shape[0]):
#     print(i)
        coord = protpdb_coord_reader(input_dat.loc[i,'Prot_PDB'], input_dat, i)
#         print('Processing {} glycan'.format(input_dat_temp.loc[i,'glycan']))
        gly = pd.read_csv(input_dat.loc[i,'glycan_pdb'], delim_whitespace=True, header=None, names=['atm','index','name','res','resid','x','y','z','f1','f2'])
        #add chain if needed
        means = gly.mean()
        gly['x'] = gly['x'] - means['x']
        gly['y'] = gly['y'] - means['y']
        gly['z'] = gly['z'] - means['z']

        dx = (float(coord.loc[0,'x2']) - float(coord.loc[0,'x1']))/3.20*10
        dy = (float(coord.loc[0,'y2']) - float(coord.loc[0,'y1']))/3.20*10
        dz = (float(coord.loc[0,'z2']) - float(coord.loc[0,'z1']))/3.20*10

        gly['x'] = gly['x'] + float(coord.iloc[0,4]) + dx*dist 
        gly['y'] = gly['y'] + float(coord.iloc[0,5]) + dy*dist
        gly['z'] = gly['z'] + float(coord.iloc[0,6]) + dz*dist

        gly = gly[0:]
        gly_string = gly.to_string(header=False,index=False,index_names=False).split('\n')
        
        with open(input_dat.loc[i,'Prot_PDB'][:-4]+'_gly.pdb', 'a') as outfile:     
            for line in (gly_string):
                pdb_line = line.split()
#             print(pdb_line)
                outfile.write("{:>4}{:>7}{:>5}{:>5}{:>5}{:>12.5}{:>8.5}{:>8.5}{:>6}{:>6}\n".format(*pdb_line))
        print('Attaching {} to N{} {}\n'.format(input_dat_temp.loc[i,'glycan']+'.pdb',input_dat.loc[i,'ASN_site'], input_dat.loc[i,'Prot_PDB'][:-4]+'_gly.pdb'))

######################## COMBINES the glycan bonds, constraints, angles, dihedrals sections to protein itps ####
def itps_combine(unique_chains, protein_path):
    for i in unique_chains:
        if os.path.isfile(protein_path[:-4]+'_'+ i +'_gly.itp'):
            os.remove(protein_path[:-4]+'_'+ i +'_gly.itp')
    
        with open (protein_path[:-4]+'_'+ i +'.itp', 'r') as infile, open (protein_path[:-4]+'_'+ i +'_gly.itp', 'a') as outfile:

            lines = (line.rstrip() for line in infile) # All lines including the blank ones
            lines = (line for line in lines if line) # Non-blank lines

            for line in lines:
                itp_line = line.split()

                if  len(itp_line) > 1 and itp_line[1] == 'bonds':
                    with open (fileDir+'/'+'gly_atoms_'+i+'.itp','r') as gly_data:
                        for line2 in gly_data:
                            outfile.write(line2)
                    outfile.write(line + '\n')

                elif len(itp_line) > 1 and itp_line[1] == 'constraints':
                    with open (fileDir+'/'+'gly_bonds_'+i+'.itp','r') as gly_data:
                        for line2 in gly_data:
                            outfile.write(line2)
                    outfile.write(line + '\n')                        

                elif len(itp_line) > 1 and itp_line[1] == 'angles':
                    if os.path.isfile('gly_constraints_'+i+'.itp'):
                        with open (fileDir+'/'+'gly_constraints_'+i+'.itp','r') as gly_data:
                            for line2 in gly_data:
                                outfile.write(line2)
                    outfile.write(line + '\n')                            

                elif len(itp_line) > 1 and itp_line[1] == 'dihedrals':
                    with open (fileDir+'/'+'gly_angles_'+i+'.itp','r') as gly_data:
                        for line2 in gly_data:
                            outfile.write(line2)
                    outfile.write(line + '\n')                        

                elif itp_line[0] == '#ifdef' and itp_line[1]=='POSRES':
                    with open (fileDir+'/'+'gly_dihedrals_'+i+'.itp','r') as gly_data:
                        for line2 in gly_data:
                            outfile.write(line2)
                    outfile.write(line + '\n')

                else:
                    outfile.write(line+'\n')

        print("Generated {}".format(protein[:-4]+'_'+ i +'_gly.itp'))
        print("Generated {}\n".format(protein[:-4]+'_'+ i +'_gly.pdb'))

    for f in glob.glob("gly_*.itp"):
        os.remove(f)
        
#######################################################################################
####################### MAIN FUNCTION ###############################################

print('\n******************* MARTINI GLYCOSYLATOR *********************\n')

unique_chains = input_dat.chain.unique()
chain_info = pd.DataFrame(dict(chain=[], atoms=[], correction=[]), dtype = int)
input_dat, chain_info = chain_splitter(unique_chains, chain_info, protein_path, protein, input_dat)

####### get indexes of asparagines 
input_dat = protpdb_reader(input_dat) 
input_dat = set_bb(input_dat, chain_info) #renumbering indexes with respect to  the chain
print("Done\n")
last_index = input_dat.loc[0,'last_index']
last_resid = input_dat.loc[0,'last_resid']

#print(input_dat)
##### GENRATING ITP FILES ##################
for i in range (0, input_dat.shape[0]):
    
#     print('processing {} {}'.format(input_dat.loc[i,'chain'],input_dat.loc[i,'Prot_PDB']))
    
    asn_par = {'bond1':[int(input_dat.loc[i,'BB']+1), last_index+2, 1, 0.478, 30000], 
            'angle1':[int(input_dat.loc[i,'BB']), int(input_dat.loc[i,'BB']+1), last_index+2, 1, 172, 170],
            'angle2':[last_index+1, last_index+2, int(input_dat.loc[i,'BB']+1), 1, 45, 50],
            'angle3':[int(input_dat.loc[i,'BB']+1), last_index+2, last_index+3, 1, 35, 120],
            'dihedral1':[int(input_dat.loc[i,'BB']), int(input_dat.loc[i,'BB']+1), last_index+2, last_index+1, 1, 80, 8, 1],
            'dihedral2':[int(input_dat.loc[i,'BB']), int(input_dat.loc[i,'BB']+1), last_index+2, last_index+3, 1, 300, 8, 1]}
    
    itp_read = itp_reader(input_dat.loc[i,'glycan_itp'])
    print("reading {}".format(input_dat_temp.loc[i,'glycan']+'.itp'))
    last_index, last_resid = itp_indexer(itp_read, last_index, last_resid, asn_par, input_dat.loc[i,'chain']) 
    print("Now indexing {} for glycan at ASN {} chain {}\n".format(input_dat_temp.loc[i,'glycan']+'.itp', input_dat.loc[i,'ASN_site'], input_dat.loc[i,'chain']))
#    Changing the last indexes and last resid of proteins based on the chain (e.g, A/B/C)
    if (i+1 != input_dat.shape[0]):
        if (input_dat.loc[i+1,'chain'] != input_dat.loc[i,'chain']):
            last_index = input_dat.loc[i+1,'last_index']
            last_resid = input_dat.loc[i+1,'last_resid']
########## ATTACHING GLYCAN PDBs TO ASN #################################################
gly_to_prot_attch(input_dat, dist, input_dat_temp)

input_dat_temp['Prot_PDB'] = input_dat['Prot_PDB']
input_dat_temp['BB_index'] = input_dat['BB']
input_dat_temp['prot_residues'] = input_dat['last_resid']
input_dat_temp['prot_atoms'] = input_dat['last_index']

print('Data used\n')
#print(input_dat.to_string(index=False)+'\n')
print(input_dat_temp.to_string(index=False)+'\n')
################ appeding glycan itp sections to protein itps ####################
itps_combine(unique_chains, protein_path)
print("\n**Successfully generated glycosylated protein pdb and topology files**\n")

################################### END #######################################

