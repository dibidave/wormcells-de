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

import datetime
dt_started = datetime.datetime.utcnow()

url = (sys.argv[1])
print(datetime.datetime.now())
print('	### ### ###  starting to process...')
print(url)

AWS_S3_ACCESS_KEY = (sys.argv[2])

AWS_S3_SECRET = (sys.argv[3])

sendgrid_key = (sys.argv[4])

sendgrid_name = (sys.argv[5])

save_path = "/home/ubuntu/"

sys.stdout = open('log_file.txt', 'w')
sys.stderr = sys.stdout
# if True:
try:
    vae_file_name = 'cpt_vae_v2.pkl'
    adata = anndata.read('cao2017packer2019taylor2019.h5ad')
    adata.var.index = adata.var['gene_id']
    adata.var.index = adata.var.index.rename('index')

    gene_dataset = GeneExpressionDataset()

    # we provide the `batch_indices` so that scvi can perform batch correction
    gene_dataset.populate_from_data(
        adata.X,
        gene_names=adata.var.index.values,
        cell_types=adata.obs['cell_type'].values,
        batch_indices=adata.obs['experiment'].cat.codes.values,
    )

    vae = VAE(gene_dataset.nb_genes, n_batch=gene_dataset.n_batches)

    trainer = UnsupervisedTrainer(
        vae,
        gene_dataset,
        train_size=0.75,
        use_cuda=False,
        frequency=1,
    )

    trainer.model.load_state_dict(torch.load(save_path + vae_file_name, map_location=torch.device('cpu')))
    print('	### ### ###  loaded vae')
    print(datetime.datetime.now())

    full = trainer.create_posterior(trainer.model, gene_dataset, indices=np.arange(len(gene_dataset)))
    latent, batch_indices, labels = full.sequential().get_latent()
    batch_indices = batch_indices.ravel()

    print('	### ### ###  computed full posterior')

    # construct the filename for saving the results
    filename = url.split('https://scvi-differential-expression.s3.us-east-2.amazonaws.com/submissions/')[1]
    filename = filename.replace('%40', '@')
    filename = filename.replace('%25', '@')
    filename = filename.replace('.csv', '')

    print('	### ### ###  url = ', url)
    # read submission csv and fetch selected cells
    submission = pd.read_csv(io.StringIO(requests.get(url).content.decode('utf-8')), index_col=0)
    selected_cells_csv_string = submission.to_csv(index=False).replace('\n', '<br>')

    # reconstruct user email from submission url
    email = url.split('https://scvi-differential-expression.s3.us-east-2.amazonaws.com/submissions/')[1]
    email = email.split('%25')[0]
    email = email.replace('%40', '@')

    # create masks for cell selection according to submission
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

    print('	### ### ###  computed idxs')

    n_samples = 5000
    M_permutation = 5000

    print('	### ### ###  computing DE...')
    print(datetime.datetime.now())
    de_res = full.differential_expression_score(
        cell_idx1.values,
        cell_idx2.values,
        n_samples=n_samples,
        mode='change',  # set to the new change mode
        M_permutation=M_permutation)

    print('	### ### ###  finished DE!')
    print(datetime.datetime.now())

    # manipulate the DE results for plotting
    de = de_res.copy()

    # we compute the ratio of the scVI scales to use that as a rough proxy for fold change
    de['ratio_scale12'] = de['scale1'] / de['scale2']
    de['log2_fold_change'] = np.log2(de['ratio_scale12'])
    # we use -log10 p-value in the volcano plot
    de['log10_pvalue'] = np.log10(de['proba_not_de'])
    # we provide the bayes factor in the CSV and on the plot mouseover
    de['abs_bayes_factor'] = np.abs(de['bayes_factor'])
    de = de.join(adata.var, how='inner')

    de['gene_color'] = 'rgba(100, 100, 100, 0.2)'
    for gene in submission['selected_genes'].dropna().values:
        gene = gene.strip()
        de['gene_color'][de['gene_name'].str.contains(gene)] = 'rgba(0, 0,255, 1)'
        de['gene_color'][de['gene_id'].str.contains(gene)] = 'rgba(0, 0,255, 1)'

    de_result_csv = de[
        ['proba_not_de', 'log10_pvalue', 'bayes_factor', 'log2_fold_change', 'gene_id','gene_name', 'gene_description']]

    de['gene_description_html'] = de['gene_description'].str.replace('\. ', '.<br>')
    print('	### ### ###  Creating plot')

    try:
        jobname = submission['job_name'][0]
    except:
        jobname = ' '

    de['gene_description_html'] = de['gene_description'].str.replace('\. ', '.<br>')
    string_bf_list = [str(bf) for bf in np.round(de['bayes_factor'].values, 3)]
    de['bayes_factor_string'] = string_bf_list

    fig = go.Figure(
        data=go.Scatter(
            x=de["log2_fold_change"].round(3)
            , y=-de["log10_pvalue"].round(3)
            , mode='markers'
            , marker=dict(color=de['gene_color'])
            , hoverinfo='text'
            , text=de['gene_description_html']
            , customdata=de.gene_id.values + '<br>Name: ' + de.gene_name.astype(str) + '<br> Bayes Factor: ' + de.bayes_factor_string
            , hovertemplate='%{customdata} <br>' +
                            '-log10(p-value): %{y}<br>' +
                            'log2 fold change: %{x}' +
                            '<extra>%{text}</extra>'
        )
        , layout={
            "title": {"text":
                          "Differential expression on C. elegans single cell data  <br> " + jobname + " <br> <a href=" + url + ">" + url + '</a>'
                , 'x': 0.5
                      }
            , 'xaxis': {'title': {"text": "log2 fold change"}}
            , 'yaxis': {'title': {"text": "-log10(p-value)"}}
        }
    )

    csv_buffer = StringIO()
    de_result_csv.to_csv(csv_buffer)
    csvfilename = 'csv/' + filename + '-results.csv'
    htmlfilename = 'plots/' + filename + '-results.html'
    print('	### ### ###  Putting files in s3...')

    client = boto3.client('s3',
                          aws_access_key_id=AWS_S3_ACCESS_KEY,
                          aws_secret_access_key=AWS_S3_SECRET
                          )

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
    print('	### ### ###  Files uploaded successfully')
    dt_ended = datetime.datetime.utcnow()
    total_time = str(int((dt_ended - dt_started).total_seconds()))

    email_body = f' Your 🌋 wormcells-de 💥 C. elegans single cell differential expression results: <br> {jobname} <br><br> <a href="{csv_url}">CSV file with results</a>  <br> <a href="{html_url}">Vocano plot</a>  <br> <a href="{url}">Your cell selection</a> <br> Processing time: {total_time}s <br> <br> <br> Thanks <br> Eduardo'
    print('	### ### ###  Email created')

    message = Mail(
        from_email='eduardo@wormbase.org',
        to_emails=email,
        subject='C. elegans single cell differential expression results',
        html_content=email_body)
    sg = SendGridAPIClient(sendgrid_key)
    response = sg.send(message)
    print(response.status_code)
    print(response.body)
    print(response.headers)
    print('	********************* Email sent ********************')
    print('	****************** Total time (s):', end = ' ')
    print(total_time)

    print('Terminating... ')
    instance_id = requests.get("http://169.254.169.254/latest/meta-data/instance-id").text

    session = boto3.Session(region_name='us-east-2',
                            aws_access_key_id=AWS_S3_ACCESS_KEY,
                            aws_secret_access_key=AWS_S3_SECRET)
    ec2 = session.resource('ec2', region_name='us-east-2')
    ec2.instances.filter(InstanceIds=[instance_id]).terminate()
# if False:
except:
    print('	XXXXXXXXXXXXXXXX SOMETHING FAILED XXXXXXXXXXXXXXX')
    filename = url.split('https://scvi-differential-expression.s3.us-east-2.amazonaws.com/submissions/')[1]
    filename = filename.replace('%40', '@')
    filename = filename.replace('%25', '@')
    filename = filename.replace('.csv', '')

    logs = open('log_file.txt', 'r').read()

    logs_buffer = StringIO()
    logs_buffer.write(logs)

    csv_buffer = StringIO()
    logsfilename = 'logs/' + filename + '-logs.txt'
    print('	### ### ###  Putting log in s3...')

    client = boto3.client('s3',
                          aws_access_key_id=AWS_S3_ACCESS_KEY,
                          aws_secret_access_key=AWS_S3_SECRET
                          )

    client.put_object(
        Body=logs_buffer.getvalue(),
        Bucket='scvi-differential-expression',
        Key=logsfilename,
        ACL='public-read'
    )
    print('	### ### ###  Log in s3. Sending email...')

    message = Mail(
        from_email='eduardo@wormbase.org',
        to_emails='veigabeltrame@gmail.com',
        subject='SOMETHING WENT WRONG!!!!111',
        html_content=logs)

    sg = SendGridAPIClient(sendgrid_key)
    response = sg.send(message)
    print(response.status_code)
    print(response.body)
    print(response.headers)

    print('	### ### ###  Email sent. Done.')

    # print('Terminating... ')
    # instance_id = requests.get("http://169.254.169.254/latest/meta-data/instance-id").text
    #
    # session = boto3.Session(region_name='us-east-2',
    #                         aws_access_key_id=AWS_S3_ACCESS_KEY,
    #                         aws_secret_access_key=AWS_S3_SECRET)
    # ec2 = session.resource('ec2', region_name='us-east-2')
    # ec2.instances.filter(InstanceIds = [instance_id]).terminate()