import os
import numpy as np
import pandas as pd
from scvi.dataset import GeneExpressionDataset
from scvi.models import VAE
from scvi.inference import UnsupervisedTrainer
import torch
import anndata
import plotly.graph_objects as go
import io
import requests
import urllib
from io import StringIO
import boto3
import datetime
import sys
from sendgrid import SendGridAPIClient
from sendgrid.helpers.mail import Mail

url = (sys.argv[1])
print(datetime.datetime.now())
print('starting to process...')
print(url)

AWS_S3_ACCESS_KEY = (sys.argv[2])

AWS_S3_SECRET = (sys.argv[3])

sendgrid_key = (sys.argv[4])

sendgrid_name = (sys.argv[5])

save_path = "/home/ubuntu/"

sys.stdout = open('log_file.txt', 'w')
sys.stderr = sys.stdout

try:
    vae_file_name = 'cpt_vae.pkl'

    adata = anndata.read('cao2017packer2019taylor2019.h5ad')

    gene_dataset = GeneExpressionDataset()

    # we provide the `batch_indices` so that scvi can perform batch correction
    gene_dataset.populate_from_data(
        adata.X,
        gene_names=adata.var.index.values,
        cell_types=adata.obs['cell_type'].values,
        batch_indices=adata.obs['experiment'].cat.codes.values,
    )

    # for this dataset 5 epochs is sufficient
    n_epochs = 10
    lr = 1e-3
    use_batches = True
    use_cuda = False

    vae = VAE(gene_dataset.nb_genes, n_batch=gene_dataset.n_batches * use_batches)

    trainer = UnsupervisedTrainer(
        vae,
        gene_dataset,
        train_size=0.75,
        use_cuda=use_cuda,
        frequency=1,
    )

    trainer.model.load_state_dict(torch.load(save_path + vae_file_name, map_location=torch.device('cpu')))
    print('loaded vae')
    print(trainer.model.eval())

    full = trainer.create_posterior(trainer.model, gene_dataset, indices=np.arange(len(gene_dataset)))
    latent, batch_indices, labels = full.sequential().get_latent()
    batch_indices = batch_indices.ravel()

    print('computed full posterior')

    filename = url.split('https://scvi-differential-expression.s3.us-east-2.amazonaws.com/submissions/')[1]
    filename = filename.replace('%40', '@')
    filename = filename.replace('%25', '@')
    filename = filename.replace('.csv', '')

    print('url = ', url)
    submission = pd.read_csv(io.StringIO(requests.get(url).content.decode('utf-8')), index_col=0)

    selected_cells_csv_string = submission.to_csv(index=False).replace('\n', '<br>')

    email = url.split('https://scvi-differential-expression.s3.us-east-2.amazonaws.com/submissions/')[1]
    email = email.split('%25')[0]
    email = email.replace('%40', '@')

    # start with all entries false
    cell_idx1 = adata.obs['cell_type'] == '000000'

    for _, entry in submission[['cell_type1', 'experiment1']].dropna().iterrows():
        cell = entry['cell_type1'].strip()
        experiment = entry['experiment1'].strip()
        curr_boolean = (adata.obs['cell_type'] == cell) & (adata.obs['experiment'] == experiment)
        cell_idx1 = (cell_idx1 | curr_boolean)

    cell_idx2 = adata.obs['cell_type'] == '000000'

    for _, entry in submission[['cell_type2', 'experiment2']].dropna().iterrows():
        cell = entry['cell_type2'].strip()
        experiment = entry['experiment2'].strip()
        curr_boolean = (adata.obs['cell_type'] == cell) & (adata.obs['experiment'] == experiment)
        cell_idx2 = (cell_idx2 | curr_boolean)

    print('computed cell_idx1 and 2')

    n_samples = 5000
    M_permutation = 5000

    print('computing DE...')
    print(datetime.datetime.now())
    de_res = full.differential_expression_score(
        cell_idx1.values,
        cell_idx2.values,
        n_samples=n_samples,
        M_permutation=M_permutation)

    print('finished DE!')
    print(datetime.datetime.now())

    # manipulate the DE results for plotting
    de = de_res.copy()

    de['ratio_mean12'] = de['mean1'] / de['mean2']
    de['ratio_scale12'] = de['scale1'] / de['scale2']

    de['log_scale_ratio'] = np.log2(de['ratio_scale12'])
    de['log_mean_ratio'] = np.log2(de['ratio_mean12'])

    de['abs_bayes_factor'] = np.abs(de['bayes1'])

    de['ratio_norm_mean12'] = de['norm_mean1'] / de['norm_mean2']
    de['log_norm_mean_ratio'] = np.log2(de['ratio_norm_mean12'])

    # make a copy of the annotated gene metadata with gene ids all lower case to avoid problems when joining dataframes
    adata_var_uppercase = adata.var.copy()
    adata_var_uppercase.index = adata_var_uppercase.index.str.upper()

    # convert top_expression gene ids index to lowercase for joining with metadata
    de = de.join(adata_var_uppercase, how='inner')

    de['gene_color'] = 'rgba(100, 100, 100, 0.2)'
    for gene in submission['selected_genes'].dropna().values:
        gene = gene.strip()
        de['gene_color'][de['gene_name'].str.contains(gene)] = 'rgba(0, 0,255, 1)'
        de['gene_color'][de['gene_id'].str.contains(gene)] = 'rgba(0, 0,255, 1)'

    de_result_csv = de[['bayes1', 'bayes2', 'log_scale_ratio', 'gene_id', 'gene_description']]

    # first we create these variables to customize the hover text in plotly's heatmap
    # the text needs to be arranged in a matrix the same shape as the heatmap
    # for the gene descriptions text, which can be several sentences, we add a line break after each sentence
    de['gene_description_html'] = de['gene_description'].str.replace('\. ', '.<br>')

    fig = go.Figure(
        data=go.Scatter(
            x=de["log_scale_ratio"].round(3)
            , y=de["abs_bayes_factor"].round(3)
            , mode='markers'
            , marker=dict(color=de['gene_color'])
            , hoverinfo='text'
            , text=de['gene_description_html']
            , customdata=de.gene_id.values + '<br>Name: ' + de.gene_name.values.astype(str)
            , hovertemplate='%{customdata} <br>' +
                            'ln(BF): %{y}<br>' +
                            'Log scale ratio: %{x}' +
                            '<extra>%{text}</extra>'
        )
        , layout={
            "title": {"text":
                          "Differential expression on C. elegans single cell data  <br><br> Original selection: <br> <small> <code> <a href=" + url + ">" + url + '</a> </small> </code>'
                , 'x': 0.5
                      }
            , 'xaxis': {'title': {
                "text": "Log2 of scVI expression scale <b>"}}
            , 'yaxis': {'title': {"text": "Natural log of absolute value of Bayes Factor, ln(BF)"}}
        }
    )

    csv_buffer = StringIO()
    de_result_csv.to_csv(csv_buffer)
    csvfilename = 'csv/' + filename + '-results.csv'
    htmlfilename = 'plots/' + filename + '-results.html'

    client = boto3.client('s3',
                          aws_access_key_id=AWS_S3_ACCESS_KEY,
                          aws_secret_access_key=AWS_S3_SECRET
                          )

    print(AWS_S3_ACCESS_KEY)
    print(AWS_S3_SECRET)
    client.put_object(
        Body=csv_buffer.getvalue(),
        Bucket='scvi-differential-expression',
        Key=csvfilename,
        ACL='public-read'
    )

    html_buffer = StringIO()
    fig.write_html(html_buffer, auto_open=True)

    client.put_object(
        Body=html_buffer.getvalue(),
        Bucket='scvi-differential-expression',
        Key=htmlfilename,
        ACL='public-read'
    )

    csv_url = 'https://scvi-differential-expression.s3.us-east-2.amazonaws.com/' + urllib.parse.quote(csvfilename)
    html_url = 'https://scvi-differential-expression.s3.us-east-2.amazonaws.com/' + urllib.parse.quote(htmlfilename)

    email_body = f' Your C. elegans single cell differential expression results from  🌋 wormcells-de 💥 are ready to download. <br><br> <a href="{csv_url}">CSV file with results</a>  <br> <a href="{html_url}">Vocano plot</a>  <br> <a href="{url}">Your original selection </a> <br> <br> Thanks <br> Eduardo'
    print(email_body)

    message = Mail(
        from_email='eduardo@wormbase.org',
        to_emails=email,
        subject='C. elegans single cell differential expression results',
        html_content=email_body)
    try:
        sg = SendGridAPIClient(sendgrid_key)
        response = sg.send(message)
        print(response.status_code)
        print(response.body)
        print(response.headers)
    except Exception as e:
        print(e.message)

    print('DONE!!!!!!!')
    print('DONE!!!!!!!')


    print('Terminating... ')
    instance_id = requests.get("http://169.254.169.254/latest/meta-data/instance-id").text

    session = boto3.Session(region_name='us-east-2',
                            aws_access_key_id=AWS_S3_ACCESS_KEY,
                            aws_secret_access_key=AWS_S3_SECRET)
    ec2 = session.resource('ec2', region_name='us-east-2')
    ec2.instances.filter(InstanceIds = [instance_id]).terminate()

except:
    str = open('log_file.txt', 'r').read()

    message = Mail(
        from_email='eduardo@wormbase.org',
        to_emails='veigabeltrame@gmail.com',
        subject='SOMETHING WENT WRONG!!!!111 user: ' + email,
        html_content=str)
    try:
        sg = SendGridAPIClient(sendgrid_key)
        response = sg.send(message)
        print(response.status_code)
        print(response.body)
        print(response.headers)
    except Exception as e:
        print(e)

    print('DONE!!!!!!!')
    print('DONE!!!!!!!')


    print('Terminating... ')
    instance_id = requests.get("http://169.254.169.254/latest/meta-data/instance-id").text

    session = boto3.Session(region_name='us-east-2',
                            aws_access_key_id=AWS_S3_ACCESS_KEY,
                            aws_secret_access_key=AWS_S3_SECRET)
    ec2 = session.resource('ec2', region_name='us-east-2')
    ec2.instances.filter(InstanceIds = [instance_id]).terminate()