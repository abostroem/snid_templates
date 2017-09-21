from astropy.io import ascii, fits
import os
import urllib.request

def parse_wiserep_table(filename):
    '''
    Copy and paste a table from WiseRep into a text file. This code parese that
    table into an astropy table so it can be searched, curated, and files can be 
    downloaded.
    '''
    names = ['Id', 'Obj. Name','Type', 'Spec Program', 'Instrument', 
             'Obs. Date', 'Observer','Red. Date','Reducer', 'ExpTime', 
             'Slit', 'Public', 'Ascii File', 'Fits File',
             'Spec Type', 'Quality', 'Remarks', 'Created By',
             'Last Modified', 'Modified By', 'junk']
    ofile = open(filename, 'r')
    all_lines = ofile.readlines()
    new_lines = []
    remove_counter=0
    for iline in all_lines:
        if "Show Files" in iline:
            all_lines.remove(iline)
    for iline in all_lines:
        if len(iline.split('\t')) < 3:
            all_lines.remove(iline)
    for even_iline, odd_iline in zip(all_lines[0::2], all_lines[1::2]):
        if odd_iline.startswith('  \t'):
            odd_iline = '\t{}'.format(odd_iline)
        even_iline = even_iline.replace(',', '-')
        odd_iline = odd_iline.replace(',', '-')
        new_lines.append(even_iline.replace(('\n'),'\t  \t').replace('\t  \t', ', ')+                          
                         odd_iline.replace('\t\t\t\n', '\n').replace('\t  \t', ', '))
        #new_lines.append(even_iline.replace(('\n'),'\t  \t')+ \
        #                 odd_iline.replace('\t\t\t\n', '\n'))
    tbdata = ascii.read(new_lines, delimiter=',', names=names, exclude_names=['junk'])
    return tbdata

def download_wiserep(tbdata, output_dir, verbose=False):
    '''
    Download files from wiserep given a table made from their search output table. Use 
    parse_wiserep_table to create this table. This code downloads fits files if they exist, 
    if not, it looks for an ascii file.
    
    If a file can't be found, it is possible I'm not searching all of the directories. Look
    at the URL, right after upload there should be a directory, if its not in wiserep_dirs, add it
    to the end or the list.
    
    TODO: This code currently only downloads final spectra. This could be an option and there 
    could be more options of other things to select on.
    '''
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    ofile = open(os.path.join(output_dir,'download_log.txt'), 'w')
    for iline in tbdata:
        if verbose:
            print(iline['Id'])
        if (iline['Quality'] == 'Final') or (iline['Quality'] == 'final'):
            base_url = 'https://wiserep.weizmann.ac.il/sites/weizmann.ac.il.wiserep/files/uploads'
            wiserep_dirs = ['/', '/', '/PESSTO/','/PESSTO-DR1/','/TS3/', '/PESSTO-DR2/']
            for idir in wiserep_dirs:
                try:
                    url = base_url+idir+'Fits/{}'.format(iline['Fits File'])
                    response = urllib.request.urlopen(url)
                    iformat='fits'
                    filename = iline['Fits File']
                    break
                except:
                    try:
                        url = base_url+idir+'Ascii/{}'.format(iline['Ascii File'])
                        response = urllib.request.urlopen(url)
                        iformat='ascii'
                        filename = iline['Ascii File']
                        break
                    except:
                        iformat = None
            if iformat is None:
                if not hasattr(iline['Fits File'], '_mask'):
                    print('Unable to locate URL for {}, ID = {}'.format(iline['Fits File'], iline['Id']))
                    ofile.write('Unable to locate URL for {}, ID = {}'.format(iline['Fits File'], iline['Id']))
                elif not hasattr(iline['Ascii File'], '_mask'):
                    print('Unable to locate URL for {}, ID = {}'.format(iline['Ascii File'], iline['Id']))
                    ofile.write('Unable to locate URL for {}, ID = {}'.format(iline['Ascii File'], iline['Id']))
                else:
                    print('No file found for ID {}'.format(iline['Id']))
                    ofile.write('No file found for ID {}\n'.format(iline['Id']))
            else:
                content = response.read()
                if iformat is 'fits':
                    f = open(os.path.join(output_dir, iline['Fits File']), "wb")
                    f.write(content)
                elif iformat is 'ascii':
                    f = open(os.path.join(output_dir, iline['Ascii File']), "w")
                    f.write(content.decode('utf-8'))
                if iformat is not None:
                    f.close()
    ofile.close()

if __name__ == "__main__":

    RAW_CAT_DIR = 'raw_catalogs'
    SAVE_CAT_DIR = 'astropy_catalogs'
    
    tbdata = parse_wiserep_table(os.path.join(RAW_CAT_DIR, '2013ej_cat_all.txt'))
    tbdata.write(os.path.join(SAVE_CAT_DIR, '2013ej_cat.csv'),  overwrite=True)
    download_wiserep(tbdata, '2013ej')
    
    tbdata = parse_wiserep_table(os.path.join(RAW_CAT_DIR, '2013by_cat_all.txt'))
    tbdata.write(os.path.join(SAVE_CAT_DIR, '2013by_cat.csv'),  overwrite=True)
    #download_wiserep(tbdata, '2013by')


    tbdata = parse_wiserep_table(os.path.join(RAW_CAT_DIR, '2012A_cat_all.txt'))
    tbdata.write(os.path.join(SAVE_CAT_DIR, '2012A_cat.csv'),  overwrite=True)
    #download_wiserep(tbdata, '2012A')


    tbdata = parse_wiserep_table(os.path.join(RAW_CAT_DIR, '2016zb_cat_all.txt'))
    tbdata.write(os.path.join(SAVE_CAT_DIR, '2016zb_cat.csv'),  overwrite=True)
    #download_wiserep(tbdata, '2016zb')






