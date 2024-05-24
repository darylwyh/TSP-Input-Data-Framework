import logging
import torch
import pickle
from torch.utils.data import Dataset
import os 
# Set up logging
logging.basicConfig(level=logging.INFO)

class TSPDataset(Dataset):
    
    def __init__(self, filename=None, size=50, num_samples=1000000, offset=0, distribution=None):
        super(TSPDataset, self).__init__()

        self.data_set = []
        if filename is not None:
            assert os.path.splitext(filename)[1] == '.pkl'

            with open(filename, 'rb') as f:
                data = pickle.load(f)
                self.data = [torch.FloatTensor(row) for row in (data[offset:offset+num_samples])]
        else:
            # Sample points randomly in [0, 1] square
            self.data = [torch.FloatTensor(size, 2).uniform_(0, 1) for i in range(num_samples)]

        self.size = len(self.data)

    def __len__(self):
        return self.size

    def __getitem__(self, idx):
        return self.data[idx]
    

# Create dataset
dataset = TSPDataset(size=50, num_samples=1000000)

# Log summary statistics
logging.info(f"Dataset size: {len(dataset)}")
logging.info(f"First item shape: {dataset[0].shape}")
logging.info(f"First item: {dataset[0]}")


'''


size = 50
data = [torch.FloatTensor(size, 2).uniform_(0, 1) for i in range(1000000)]
size = len(data)

print(data)
print("size is: ", size)


'''