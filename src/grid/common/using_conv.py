
import os
import sys

if len(sys.argv) < 2:
    print('Must specify file to process')
    quit()

file_to_proc = sys.argv[1]
orig_file    = file_to_proc+'.orig'

os.rename(file_to_proc, orig_file)


with open(orig_file,'r') as fin:
    lines = fin.readlines()

with open(file_to_proc,'w') as fout:
    
    for line in lines:
        line_out = line
        sline = line.split();    
        if len(sline) == 0:
            fout.write('\n')
            continue
        if (sline[0] == 'typedef') and (sline[-1][-1] == ';'):

            sline.remove('typedef')
            
            typename = ''
            if sline[0] == 'typename':
                typename = ' typename'
                sline.remove('typename')

            variable = sline[-1].replace(';','')
            sline.pop()
  
            actual = ''.join(sline) 

            new_line = 'using '+variable+' ='+typename+' '+actual+';\n'
            print(' ')
            print('<< '+line_out)
            print('>> '+new_line)                
            ans = input('Replace? ')
            if ans in ['y','Y']:
                line_out = new_line

        fout.write(line_out)
