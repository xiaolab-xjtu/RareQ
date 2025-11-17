from collections import Counter
import matplotlib.pyplot as plt
from pysankey2 import Sankey
from scipy import sparse
from scipy.sparse import hstack, vstack, coo_matrix
from scipy.io import mmread
from sklearn import metrics
import math
import time
import anndata
import torch
import numpy as np
import pandas as pd
import scanpy as sc
import random
import os
from warnings import filterwarnings
import scipy.sparse as sp
from operator import itemgetter
filterwarnings("ignore")
import json, os
import math, copy, time
import numpy as np
from collections import defaultdict
import pyHGT2
from pyHGT2.data import *
from pyHGT2.utils import *
from pyHGT2.model1 import *
from pyHGT2.conv import *
from sklearn import metrics

from collections import Counter
import matplotlib.pyplot as plt
from pysankey2 import Sankey
from scipy.sparse import hstack, vstack, coo_matrix
from scipy.io import mmread
from sklearn import metrics
import math
import time
import anndata
import torch
import numpy as np
import pandas as pd
import scanpy as sc
import random
import os
from scipy import sparse
from scipy.io import mmread
from scipy.sparse import hstack, vstack, coo_matrix

seed = 0
torch.manual_seed(seed)
torch.cuda.manual_seed(seed)
torch.cuda.manual_seed_all(seed)
random.seed(seed)
np.random.seed(seed)
os.environ['PYTHONHASHSEED'] = str(seed)
torch.manual_seed(seed)
torch.cuda.manual_seed(seed)


## Sim-PBMC 1-3

os.chdir('./Simulation/')

pbmc = sc.read_h5ad('GSE194122_openproblems_neurips2021_multiome_BMMC_processed.h5ad')
cell_type = pbmc.obs

## Gene_Peak.mtx, Gene_Cell.mtx, Peak_Cell.mtx, Gene_names.tsv, Cell_names.tsv, Peak_names.tsv are downloaded from zenodo: https://zenodo.org/records/8406470

gene_peak = anndata.read_mtx('Gene_Peak.mtx')
gene_cell = anndata.read_mtx('Gene_Cell.mtx')
peak_cell = anndata.read_mtx('Peak_Cell.mtx')
gene_names = pd.read_csv('Gene_names.tsv', sep='\t', header=None)
cell_names = pd.read_csv('Cell_names.tsv', sep='\t', header=None)
peak_names = pd.read_csv('Peak_names.tsv', sep='\t', header=None)
peak_cell.obs_names = peak_names[0]
peak_cell.var_names = cell_names[0]
gene_cell.obs_names = gene_names[0]
gene_cell.var_names = cell_names[0]
gene_peak.obs_names = gene_names[0]
gene_peak.var_names = peak_names[0]

#Gene_Cell, Peak_Cell, Gene_Peak = gene_cell,peak_cell,gene_peak
true_label = cell_type.loc[cell_names[0].tolist(), :]
cell_type['number'] = 1
cell_count = cell_type.groupby('cell_type').count()
cell_count =cell_count.sort_values(by='number',ascending=True, inplace=False)

true_label.to_csv('site1/True_label_site1.csv')

cell_count['percentage'] = (cell_count['number'] / cell_count['number'].sum())
cell_count

all(true_label.index == cell_names[0])

# Reset index
cell_df = true_label.reset_index(drop=True)

# Group by cell type to get the corresponding index list
cell_type_dict = cell_df.groupby('cell_type').apply(lambda x: x.index.tolist()).to_dict()

# Sort the dictionary based on the number of cells in each type
sorted_cell_type_dict = {k: v for k, v in sorted(cell_type_dict.items(), key=lambda item: len(item[1]))}

# Get the lengths of each cell type's indices
lengths = [len(indices) for indices in sorted_cell_type_dict.values()]

# Calculate the total count
total = sum(lengths)

# Print the proportion of each cell type
for cell_type, indices in sorted_cell_type_dict.items():
    proportion = len(indices) / total
    print(f"Cell type: {cell_type}, Proportion: {proportion:.3f}")
    
# Uncomment the lines below to print the keys and lengths of the sorted dictionary
# print(sorted_cell_type_dict.keys())
# for cell_type, indices in sorted_cell_type_dict.items():
#     print(f"Cell type: {cell_type}, Length: {len(indices)}")


# Sim-PBMC 1
import numpy as np
import random

def sample_cells(cell_type_dict, total_samples=5000):
    # Extract all cell types
    cell_types = list(cell_type_dict.keys())
    
    # The last 4 types must be sampled
    fixed_types = cell_types[-4:]
    
    # Randomly select 1 cell type from the first types
    random_type = random.choice(cell_types[:2])
    
    # Combine these 5 types
    sample_types = fixed_types + [random_type]
    
    # Calculate the proportion of each cell type within the selected types
    type_counts = {cell_type: len(cells) for cell_type, cells in cell_type_dict.items() if cell_type in sample_types}
    total_count = sum(type_counts.values()) - 50 if random_type in type_counts else 0
    type_ratio = {cell_type: count / total_count for cell_type, count in type_counts.items() if cell_type != random_type}
    
    # Sample these cell types proportionally to a total of 5000
    sample_counts = {cell_type: int(np.round(ratio * total_samples)) for cell_type, ratio in type_ratio.items()}
    
    # Find the difference between the total number of samples and 5000 (including the fixed 50), then allocate this difference to the cell type with the largest sample count
    diff = total_samples - sum(sample_counts.values()) - 50
    max_type = max(sample_counts, key=sample_counts.get)
    sample_counts[max_type] += diff
    
    # Perform sampling
    sample_indices = []
    for cell_type, count in sample_counts.items():
        cells = cell_type_dict[cell_type]
        if len(cells) > count:
            cells = np.random.choice(cells, size=count, replace=False)
        sample_indices.extend(cells)
    
    # Additionally sample 50 samples from the selected random type
    cells = cell_type_dict[random_type]
    if len(cells) > 50:
        cells = np.random.choice(cells, size=50, replace=False)
    sample_indices.extend(cells)
        
    return sample_indices

import os
import numpy as np
import scipy.sparse
from tqdm import tqdm
from scipy import io


ip_file = '/data/Home/fabotao/Projects/Rare_cell/data/PBMCs/site1/'
filetype_name_list = ['Rare1_Ordinary4']

for filetype_name in filetype_name_list:
    for file_num in tqdm(range(0, 50)):  # Process the first 50 folders
        path2 = os.path.join(ip_file, filetype_name, str(file_num), 'R_input')
        
        # Set random seed
        np.random.seed(file_num)
        random.seed(file_num)
        
        # Execute sampling function
        sample_indices = sample_cells(sorted_cell_type_dict, total_samples=5000)
        
        # Sample corresponding gene_cell, peak_cell matrices, and true_label labels
        sub_gene_cell = gene_cell[:, sample_indices]
        sub_peak_cell = peak_cell[:, sample_indices]
        gene_cell_matrix = sub_gene_cell.X
        peak_cell_matrix = sub_peak_cell.X
        
        sub_label = list(true_label.iloc[sample_indices]['cell_type'])
        
        print('sub_gene_cell.shape:', end=' ')
        print(sub_gene_cell.shape)
        print('sub_peak_cell.shape:', end=' ')
        print(sub_peak_cell.shape)
        
        # Save files
        io.mmwrite(os.path.join(path2, "sub_gene_cell.mtx"), gene_cell_matrix)
        io.mmwrite(os.path.join(path2, "sub_peak_cell.mtx"), peak_cell_matrix)
        with open(os.path.join(path2, "sub_label.txt"), "w") as file:
            for label in sub_label:
                file.write(f"{label}\n")
        with open(os.path.join(path2, "sample_indices.txt"), "w") as file:
            for index in sample_indices:
                file.write(f"{index}\n")
                
                

# Sim-PBMC 2

import numpy as np
import random

def sample_cells(cell_type_dict, total_samples=5000):
    # Extract all cell types
    cell_types = list(cell_type_dict.keys())
    
    # The last 9 types must be sampled
    fixed_types = cell_types[-9:]
    
    # Randomly select 1 cell type from the first types
    random_type = random.choice(cell_types[:2])
    
    # Combine these 10 types
    sample_types = fixed_types + [random_type]
    
    # Calculate the proportion of each cell type within the selected types
    type_counts = {cell_type: len(cells) for cell_type, cells in cell_type_dict.items() if cell_type in sample_types}
    total_count = sum(type_counts.values()) - 50 if random_type in type_counts else 0
    type_ratio = {cell_type: count / total_count for cell_type, count in type_counts.items() if cell_type != random_type}
    
    # Sample these cell types proportionally to a total of 5000
    sample_counts = {cell_type: int(np.round(ratio * total_samples)) for cell_type, ratio in type_ratio.items()}
    
    # Find the difference between the total number of samples and 5000 (including the fixed 50), then allocate this difference to the cell type with the largest sample count
    diff = total_samples - sum(sample_counts.values()) - 50
    max_type = max(sample_counts, key=sample_counts.get)
    sample_counts[max_type] += diff
    
    # Perform sampling
    sample_indices = []
    for cell_type, count in sample_counts.items():
        cells = cell_type_dict[cell_type]
        if len(cells) > count:
            cells = np.random.choice(cells, size=count, replace=False)
        sample_indices.extend(cells)
    
    # Additionally sample 50 samples from the selected random type
    cells = cell_type_dict[random_type]
    if len(cells) > 50:
        cells = np.random.choice(cells, size=50, replace=False)
    sample_indices.extend(cells)
        
    return sample_indices

import os
import numpy as np
import scipy.sparse
from tqdm import tqdm
from scipy import io

ip_file = '/data/Home/fabotao/Projects/Rare_cell/data/PBMCs/site1/'
filetype_name_list = ['Rare1_Ordinary9']

for filetype_name in filetype_name_list:
    for file_num in tqdm(range(0, 50)):  # Process the first 50 folders
        path2 = os.path.join(ip_file, filetype_name, str(file_num), 'R_input')
        
        # os.makedirs(path2)
        # Set random seed
        np.random.seed(file_num)
        random.seed(file_num)
        
        # Execute sampling function
        sample_indices = sample_cells(sorted_cell_type_dict, total_samples=5000)
        
        # Sample corresponding gene_cell, peak_cell matrices, and true_label labels
        sub_gene_cell = gene_cell[:, sample_indices]
        sub_peak_cell = peak_cell[:, sample_indices]
        gene_cell_matrix = sub_gene_cell.X
        peak_cell_matrix = sub_peak_cell.X
        
        sub_label = list(true_label.iloc[sample_indices]['cell_type'])
        
        print('sub_gene_cell.shape:', end=' ')
        print(sub_gene_cell.shape)
        print('sub_peak_cell.shape:', end=' ')
        print(sub_peak_cell.shape)
        
        # Save files
        io.mmwrite(os.path.join(path2, "sub_gene_cell.mtx"), gene_cell_matrix)
        io.mmwrite(os.path.join(path2, "sub_peak_cell.mtx"), peak_cell_matrix)
        with open(os.path.join(path2, "sub_label.txt"), "w") as file:
            for label in sub_label:
                file.write(f"{label}\n")
        with open(os.path.join(path2, "sample_indices.txt"), "w") as file:
            for index in sample_indices:
                file.write(f"{index}\n")
                



# Sim-PBMC 3
import numpy as np
import random

def sample_cells(cell_type_dict, total_samples=5000):
    # Extract all cell types
    cell_types = list(cell_type_dict.keys())
    
    # The last ten types must be sampled
    fixed_types = cell_types[-10:]
    
    # Randomly select 5 cell types from the first types
    random_types = random.sample(cell_types[:8], 5)
    
    # Combine these 15 types
    sample_types = fixed_types + random_types
    
    # Calculate the proportion of each cell type within these 15 types
    type_counts = {cell_type: len(cells) for cell_type, cells in cell_type_dict.items() if cell_type in sample_types}
    total_count = sum(type_counts.values())
    type_ratio = {cell_type: count / total_count for cell_type, count in type_counts.items()}
    
    # Sample these cell types proportionally to a total of 5000
    sample_counts = {cell_type: int(np.round(ratio * total_samples)) for cell_type, ratio in type_ratio.items()}
    
    # Find the difference between the total number of samples and 5000, then allocate this difference to the cell type with the largest sample count
    diff = total_samples - sum(sample_counts.values())
    max_type = max(sample_counts, key=sample_counts.get)
    sample_counts[max_type] += diff
    
    # Perform sampling
    sample_indices = []
    for cell_type, count in sample_counts.items():
        cells = cell_type_dict[cell_type]
        if len(cells) > count:
            cells = np.random.choice(cells, size=count, replace=False)
        sample_indices.extend(cells)
        
    return sample_indices
  
import os
import numpy as np
import scipy.sparse
from tqdm import tqdm
from scipy import io

ip_file = '/data/Home/fabotao/Projects/Rare_cell/data/PBMCs/site1/'
filetype_name_list = ['Rare5_Ordinary10']

for filetype_name in filetype_name_list:
    for file_num in tqdm(range(0, 50)):  # Process 50 folders
        path2 = os.path.join(ip_file, filetype_name, str(file_num), 'R_input')
        # os.makedirs(path2)
        
        # Set random seed
        np.random.seed(file_num)
        random.seed(file_num)
        
        # Execute sampling function
        sample_indices = sample_cells(sorted_cell_type_dict, total_samples=5000)
        
        # Sample corresponding gene_cell, peak_cell matrices, and true_label labels
        sub_gene_cell = gene_cell[:, sample_indices]
        sub_peak_cell = peak_cell[:, sample_indices]
        gene_cell_matrix = sub_gene_cell.X
        peak_cell_matrix = sub_peak_cell.X
        
        sub_label = list(true_label.iloc[sample_indices]['cell_type'])
        
        print('sub_gene_cell.shape:', end=' ')
        print(sub_gene_cell.shape)
        print('sub_peak_cell.shape:', end=' ')
        print(sub_peak_cell.shape)
        
        # Save files
        io.mmwrite(os.path.join(path2, "sub_gene_cell.mtx"), gene_cell_matrix)
        io.mmwrite(os.path.join(path2, "sub_peak_cell.mtx"), peak_cell_matrix)
        with open(os.path.join(path2, "sub_label.txt"), "w") as file:
            for label in sub_label:
                file.write(f"{label}\n")
        with open(os.path.join(path2, "sample_indices.txt"), "w") as file:
            for index in sample_indices:
                file.write(f"{index}\n")
                

