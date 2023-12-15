import sys
import os
import pandas as pd


wdl_output_dirs = sys.argv[1:]
file_names = [f for f in os.listdir(wdl_output_dirs[0]) if (f != 'README.txt') and (f[0] != '.')]

# Make output dir location
print("Creating directory for WDL outputs...")
try:
    os.mkdir('./wdl_outputs')
except FileExistsError:
    print("ERROR: The directory wdl_outputs already exists. Please delete or archive it before rerunning.")
    sys.exit(1)

for file in file_names:
    print(f'Combining file {file} across different directories...')
    full_df = pd.DataFrame()
    for output_dir in wdl_output_dirs:
        print(f'Loading {file} from {output_dir}...')
        df = pd.read_csv(f'{output_dir}/{file}', sep='\t', low_memory=False)
        full_df = pd.concat([full_df, df])
    full_df.to_csv(f'wdl_outputs/{file}', sep='\t', index=False)

print('Done!')