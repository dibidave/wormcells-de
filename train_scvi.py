import os
from scvi.dataset import GeneExpressionDataset
from scvi.models import VAE
from scvi.inference import UnsupervisedTrainer
import torch
import anndata
import pandas as pd

save_path = '/home/ubuntu/'
vae_file_name = 'worm_vae_scvi_061.pkl'
adata_filename='cao2017packer2019taylor2019.h5ad'
adata = anndata.read(adata_filename)

gene_dataset = GeneExpressionDataset()

# we provide the `batch_indices` so that scvi can perform batch correction
gene_dataset.populate_from_data(
            adata.X,
            gene_names=adata.var.index.values,
            cell_types=adata.obs['cell_type'].values,
            batch_indices=adata.obs['experiment'].cat.codes.values,
            )

n_epochs = 5
lr = 1e-3

vae = VAE(gene_dataset.nb_genes, n_batch=gene_dataset.n_batches)

trainer = UnsupervisedTrainer(
    vae,
    gene_dataset,
    train_size=0.75,  # number between 0 and 1, default 0.8
    frequency=1
)

full_file_save_path = os.path.join(save_path, vae_file_name)
print('Training model and saving on path:')
print(full_file_save_path)
trainer.train(n_epochs=n_epochs, lr=lr)
torch.save(trainer.model.state_dict(), full_file_save_path)


train_test_results = pd.DataFrame(trainer.history).rename(columns={'elbo_train_set':'Train', 'elbo_test_set':'Test'})

print(train_test_results)

print('******** !!! DONE !!! ******')