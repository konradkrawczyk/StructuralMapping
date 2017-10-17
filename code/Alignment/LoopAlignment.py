from FREAD.pyfread_api import run_fread
if __name__ == '__main__':

	loop = 'H2'

	run_fread('../../data/fread_db/db_CDR'+loop,'12e8.pdb',51,'DPEIGDDD','P','summary.txt')
